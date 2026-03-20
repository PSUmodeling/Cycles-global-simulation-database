import argparse
import cycles.gadm as gadm
import cycles.weather as weather
import geopandas as gpd
import numpy as np
import pandas as pd
import rioxarray
import xarray
from dataclasses import dataclass, field
from pathlib import Path
from shapely.geometry import Point
from config import CROPGRIDS, GAEZ
from config import DATA_DIR, TEMP_DIR, LUT_CSV
from config import WGS84, CEA

GAEZ_NAMES = {
    'bean': 'Pulses',
    'cassava': 'Cassava',
    'lentil': 'Pulses',
    'maize': 'Maize',
    'millet': 'Millet',
    'potato': 'PotatoAndSweetpotato',
    'rice': 'Rice',
    'sorghum': 'Sorghum',
    'soybean': 'Soybean',
    'sweetpotato': 'PotatoAndSweetpotato',
    'wheat': 'Wheat',
}

GAEZ_VARIABLES = {
    'harvested_area': 'HarvArea',
    'production': 'Production',
}

@dataclass
class Crop:
    name: str
    data: pd.DataFrame = field(init=False)

    def __post_init__(self) -> None:
        nc = xarray.open_dataset(CROPGRIDS.file_path(self.name))
        nc.rio.write_crs(WGS84, inplace=True)
        df = pd.DataFrame(nc['croparea'].to_series().rename('crop_area_ha'))
        df.replace(-1, np.nan, inplace=True)

        for management in ['Rainfed', 'Irrigated']:
            for var, gaez_var in GAEZ_VARIABLES.items():
                rds = rioxarray.open_rasterio(GAEZ.file_path(GAEZ_NAMES[self.name], management, gaez_var), masked=True)
                rds = rds.rio.reproject_match(nc)
                df = pd.concat([df, rds.to_series().rename(f'{management.lower()}_{var}').droplevel(level='band')], axis=1)

        df['rainfed_fraction'] = df['rainfed_harvested_area'] / (df['rainfed_harvested_area'] + df['irrigated_harvested_area'])
        df['rainfed_fraction'] = df['rainfed_fraction'].fillna(1.0)

        df['rainfed_crop_area_ha'] = df['crop_area_ha'] * df['rainfed_fraction']
        df['irrigated_crop_area_ha'] = df['crop_area_ha'] - df['rainfed_crop_area_ha']
        df['grid_index'] = range(len(df))
        self.data = df.reset_index(names=['latitude', 'longitude'])


def get_country_gdf(global_gdf: gpd.GeoDataFrame, country: str) -> gpd.GeoDataFrame:
    gdf = global_gdf[global_gdf['GID_0'] == country].copy().to_crs(CEA)
    gdf['region_area_km2'] = gdf.area / 1.0E6
    gdf['centroid'] = gdf.centroid

    return gdf


def get_cropgrids_gdf(crop_data: pd.DataFrame, gadm_gdf: gpd.GeoDataFrame, country: str) -> gpd.GeoDataFrame:
    df = pd.read_csv(TEMP_DIR / f'{country}.csv')

    if df.empty: return gpd.GeoDataFrame()

    # Calculate fraction weighted crop area and production for each grid cell
    for col in ['rainfed_crop_area_ha', 'irrigated_crop_area_ha', 'rainfed_production', 'irrigated_production']:
        df[col] = df[['grid_index', 'grid_fraction']].apply(lambda x: crop_data.iloc[int(x['grid_index'])][col] * x['grid_fraction'], axis=1)

    # Filter out grid cells with zero crop area to speed up distance calculations
    df = df[(df['rainfed_crop_area_ha'] > 0.0) | (df['irrigated_crop_area_ha'] > 0.0)]

    if df.empty: return gpd.GeoDataFrame()

    df['latlon'] = df['grid_index'].map(lambda x: (crop_data.iloc[int(x)]['latitude'], crop_data.iloc[int(x)]['longitude']))
    df['geometry'] = df['latlon'].map(lambda x: Point(x[1], x[0]))

    gdf = gpd.GeoDataFrame(df, geometry=df['geometry'], crs=WGS84).to_crs(CEA)
    gdf.set_index('GID', inplace=True)

    gdf = gdf.join(gadm_gdf['centroid'])
    gdf['distance_m'] = gdf['centroid'].distance(gdf['geometry'])

    return gdf


def find_reference_coordinate(gid, gadm_gdf, grid_gdf, management):
    sub_gdf = grid_gdf.loc[[gid]]
    sub_gdf = sub_gdf[sub_gdf.within(gadm_gdf.loc[gid]['geometry'])]

    if sub_gdf.empty:
        region_centroid = gadm_gdf.loc[[gid]]['centroid'].to_crs(WGS84).iloc[0]
        print(f"No grid cell with crop area found within the region {gid}. Using region centroid as reference coordinate.")
        return (region_centroid.y, region_centroid.x)
    else:
        sub_gdf = sub_gdf.sort_values(by=[f'{management}_crop_area_ha', f'{management}_production', 'distance_m'], ascending=[False, False, True])

        return sub_gdf.iloc[0]['latlon']


def area_fraction_filter(row):
    return (row['grid_fraction'] > 0.5) & ((row['crop_area_ha'] > 0.1) | (row['crop_area_fraction'] > 1.0E-3))


def main(crop_name):
    # Read CROPGRIDS and GAEZ data for selected crop
    crop = Crop(crop_name)

    # Read global administrative boundaries
    global_gdf = gadm.read_gadm(DATA_DIR / 'gadm', 'global', 'county')

    lookup_dfs = {
        'rainfed': pd.DataFrame(),
        'irrigated': pd.DataFrame(),
    }

    for gid in global_gdf['GID_0'].unique():
        gadm_gdf = get_country_gdf(global_gdf, gid)
        grid_gdf = get_cropgrids_gdf(crop.data, gadm_gdf, gid)

        if grid_gdf.empty: continue

        _dfs = {}

        for management in ['rainfed', 'irrigated']:
            _dfs[management] = grid_gdf.reset_index()[['GID', f'{management}_crop_area_ha', 'grid_fraction']].groupby('GID').agg(
                {
                    f'{management}_crop_area_ha': 'sum',
                    'grid_fraction': 'max',
                }
            )
            _dfs[management] = gadm_gdf[['UID', 'NAME_0', 'NAME_1', 'NAME_2', 'region_area_km2']].join(_dfs[management], how='inner')
            _dfs[management].rename(columns={f'{management}_crop_area_ha': 'crop_area_ha'}, inplace=True)
            _dfs[management]['crop_area_fraction'] = _dfs[management]['crop_area_ha'] * 1.0E-2 / _dfs[management][f'region_area_km2']

            _dfs[management] = _dfs[management][area_fraction_filter(_dfs[management])]

            if not _dfs[management].empty:
                _dfs[management]['reference_coordinate'] = _dfs[management].apply(lambda x: find_reference_coordinate(x.name, gadm_gdf, grid_gdf, management), axis=1)

                locations = _dfs[management]['reference_coordinate'].tolist()
                _dfs[management]['weather'] = weather.find_grids('GLDAS', locations=locations, screen_output=False, remove_duplicates=False)
                _dfs[management]['soil'] = _dfs[management].index.map(lambda x: f'{crop.name}_{management}_soilgrids_{x}.soil')

                lookup_dfs[management] = pd.concat([lookup_dfs[management], _dfs[management].drop(columns=['grid_fraction', 'crop_area_fraction'])], axis=0)

    for management in ['rainfed', 'irrigated']:
        lookup_dfs[management]['reference_latitude'] = lookup_dfs[management]['reference_coordinate'].map(lambda x: x[0])
        lookup_dfs[management]['reference_longitude'] = lookup_dfs[management]['reference_coordinate'].map(lambda x: x[1])
        lookup_dfs[management].sort_values(by=['UID'], inplace=True)
        lookup_dfs[management].drop(columns=['UID', 'reference_coordinate']).to_csv(LUT_CSV(crop.name, management, 'global'), float_format='%.3f')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Create crop lookup files for global administrative regions")
    parser.add_argument(
        '--crop',
        choices=list(GAEZ_NAMES.keys()),
        help='Crop',
    )
    args = parser.parse_args()

    print(f"Generating global crop look-up tables for {args.crop}...")
    main(args.crop)
