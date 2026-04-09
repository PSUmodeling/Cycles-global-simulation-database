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
from config import CROPS, MANAGEMENTS
from config import SOIL_SOURCE, WEATHER_SOURCE
from config import CROPGRIDS, GAEZ, GADM
from config import TEMP_DIR, LUT_CSV
from config import WGS84, CEA

AREA_FRACTION_THRESHOLD = 1e-3
GRID_FRACTION_THRESHOLD = 0.5
MIN_CROP_AREA_HA = 0.01

def _load_gaez_series(crop_name: str, management: str, reference_nc: xarray.Dataset) -> pd.DataFrame:
    """Load GAEZ harvested-area and production rasters and reproject to the CROPGRIDS grid, returning a DataFrame with one column per variable.
    """
    frames: dict[str, pd.Series] = {}
    for var in ['harvested_area', 'production']:
        rds = rioxarray.open_rasterio(GAEZ.path(crop_name, management, var), masked=True)
        series = rds.rio.reproject_match(reference_nc).to_series().rename(f'{management.lower()}_{var}').droplevel('band')
        frames[f'{management.lower()}_{var}'] = series
    return pd.DataFrame(frames)


@dataclass
class Crop:
    name: str
    data: pd.DataFrame = field(init=False)

    def __post_init__(self) -> None:
        self.data = _load_crop_data(self.name)


def _load_crop_data(crop_name: str) -> pd.DataFrame:
    """Load CROPGRIDS + GAEZ data for one crop and return a tidy DataFrame.
    """
    with xarray.open_dataset(CROPGRIDS.path(crop_name)) as nc:
        nc.rio.write_crs(WGS84, inplace=True)
        df = pd.DataFrame({
            'crop_area_ha': nc['croparea'].to_series(),
            'harvested_area_ha': nc['harvarea'].to_series(),
        })

    df.replace(CROPGRIDS.ocean_value, np.nan, inplace=True)

    for management in MANAGEMENTS:
        df = pd.concat([df, _load_gaez_series(crop_name, management, nc)], axis=1)

    rainfed_ha = df['rainfed_harvested_area']
    irrigated_ha = df['irrigated_harvested_area']
    df['rainfed_fraction'] = (rainfed_ha / (rainfed_ha + irrigated_ha)).fillna(1.0)

    df['rainfed_crop_area_ha'] = df['crop_area_ha'] * df['rainfed_fraction']
    df['irrigated_crop_area_ha'] = df['crop_area_ha'] - df['rainfed_crop_area_ha']
    df['rainfed_harvested_area_ha'] = df['harvested_area_ha'] * df['rainfed_fraction']
    df['irrigated_harvested_area_ha'] = df['harvested_area_ha'] - df['rainfed_harvested_area_ha']
    df['grid_index'] = range(len(df))

    return df.reset_index(names=['latitude', 'longitude'])


def get_country_gdf(global_gdf: gpd.GeoDataFrame, country: str) -> gpd.GeoDataFrame:
    gdf = global_gdf[global_gdf['GID_0'] == country].copy().to_crs(CEA)
    gdf['region_area_km2'] = gdf.area / 1e6
    gdf['centroid'] = gdf.centroid
    return gdf


def get_cropgrids_gdf(crop_data: pd.DataFrame, gadm_gdf: gpd.GeoDataFrame, country: str) -> gpd.GeoDataFrame:
    """Load the pre-computed country intersection CSV, scale crop columns by grid_fraction, and return a GeoDataFrame indexed by GID.
    """
    df = pd.read_csv(TEMP_DIR / f'{country}.csv')
    if df.empty:
        return gpd.GeoDataFrame()

    # Vectorized fraction-weighted scaling
    scaled_cols = [
        'rainfed_crop_area_ha', 'irrigated_crop_area_ha',
        'rainfed_harvested_area_ha', 'irrigated_harvested_area_ha',
        'rainfed_production', 'irrigated_production',
    ]
    source = crop_data.iloc[df['grid_index'].astype(int).values].reset_index(drop=True)
    for col in scaled_cols:
        df[col] = source[col].values * df['grid_fraction'].values

    df = df[(df['rainfed_crop_area_ha'] > 0) | (df['irrigated_crop_area_ha'] > 0)]
    if df.empty:
        return gpd.GeoDataFrame()

    # Build geometry directly from crop_data without an intermediate tuple column
    df['latitude'] = source.loc[df.index, 'latitude'].values
    df['longitude'] = source.loc[df.index, 'longitude'].values
    df['latlon'] = list(zip(df['latitude'], df['longitude']))
    df['geometry'] = gpd.points_from_xy(df['longitude'], df['latitude'])

    gdf = gpd.GeoDataFrame(df, geometry='geometry', crs=WGS84).to_crs(CEA)
    gdf = gdf.set_index('GID').join(gadm_gdf['centroid'])
    gdf['distance_m'] = gdf['centroid'].distance(gdf['geometry'])

    return gdf


def find_reference_coordinate(gid: str, gadm_gdf: gpd.GeoDataFrame, grid_gdf: gpd.GeoDataFrame, management: str) -> tuple[float, float]:
    """Return (lat, lon) of the best representative grid cell for a region.
    """
    sub_gdf = grid_gdf.loc[[gid]]
    sub_gdf = sub_gdf[sub_gdf.within(gadm_gdf.loc[gid, 'geometry'])]

    if sub_gdf.empty:
        # Fallback: centroid is in CEA; convert to WGS84 for lat/lon
        centroid_wgs84 = gpd.GeoSeries([gadm_gdf.loc[gid, 'centroid']], crs=CEA).to_crs(WGS84).iloc[0]
        print(f"No grid cell with crop area found within {gid}. Using region centroid as reference coordinate.")
        return (centroid_wgs84.y, centroid_wgs84.x)

    best = sub_gdf.sort_values(
        by=[f'{management}_crop_area_ha', f'{management}_production', 'distance_m'],
        ascending=[False, False, True],
    )
    return best.iloc[0]['latlon']


def _passes_area_fraction_filter(df: pd.DataFrame) -> pd.Series:
    return (df['grid_fraction'] > GRID_FRACTION_THRESHOLD) & ((df['crop_area_ha'] > MIN_CROP_AREA_HA) | (df['crop_area_fraction'] > AREA_FRACTION_THRESHOLD))


def build_management_df(management: str, gadm_gdf: gpd.GeoDataFrame, grid_gdf: gpd.GeoDataFrame, crop_name: str) -> pd.DataFrame:
    """Aggregate grid-cell data to region level for one management type and apply the area-fraction filter.
    Returns an empty DataFrame if nothing passes the filter.
    """
    area_col = f'{management}_crop_area_ha'
    harvested_col = f'{management}_harvested_area_ha'

    agg = grid_gdf.reset_index()[['GID', area_col, harvested_col, 'grid_fraction']].groupby('GID').agg({area_col: 'sum', harvested_col: 'sum', 'grid_fraction': 'max'})

    df = gadm_gdf[['UID', 'NAME_0', 'NAME_1', 'NAME_2', 'region_area_km2']].join(agg, how='inner')
    df.rename(columns={area_col: 'crop_area_ha', harvested_col: 'harvested_area_ha'}, inplace=True)
    df['crop_area_fraction'] = df['crop_area_ha'] * 1e-2 / df['region_area_km2']

    df = df[_passes_area_fraction_filter(df)]
    if df.empty:
        return df

    df['reference_coordinate'] = df.apply(lambda row: find_reference_coordinate(row.name, gadm_gdf, grid_gdf, management), axis=1)
    locations = df['reference_coordinate'].tolist()
    df['weather'] = weather.find_grids(WEATHER_SOURCE, locations=locations, screen_output=False, remove_duplicates=False)
    df['soil'] = df.index.map(lambda gid: f'{crop_name}_{management}_{SOIL_SOURCE}_{gid}.soil')

    return df.drop(columns=['grid_fraction', 'crop_area_fraction'])


def main(crop_name: str) -> None:
    crop = Crop(crop_name)
    global_gdf = gadm.read_gadm(GADM.path, 'global', 'county')

    lookup_dfs: dict[str, list[pd.DataFrame]] = {m: [] for m in MANAGEMENTS}

    for country in global_gdf['GID_0'].unique():
        if not (TEMP_DIR / f'{country}.csv').exists():
            continue

        gadm_gdf = get_country_gdf(global_gdf, country)
        grid_gdf = get_cropgrids_gdf(crop.data, gadm_gdf, country)

        if grid_gdf.empty:
            continue

        for management in MANAGEMENTS:
            df = build_management_df(management, gadm_gdf, grid_gdf, crop_name)
            if not df.empty:
                lookup_dfs[management].append(df)

    for management in MANAGEMENTS:
        if not lookup_dfs[management]:
            continue
        combined = pd.concat(lookup_dfs[management], axis=0)
        combined['reference_latitude']  = combined['reference_coordinate'].map(lambda x: x[0])
        combined['reference_longitude'] = combined['reference_coordinate'].map(lambda x: x[1])
        combined.sort_values('UID', inplace=True)
        LUT_CSV(crop_name, management, 'global').parent.mkdir(parents=True, exist_ok=True)
        combined.drop(columns=['UID', 'reference_coordinate']).to_csv(LUT_CSV(crop_name, management, 'global'), float_format='%.3f')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Create crop lookup files for global administrative regions")
    parser.add_argument('--crop', choices=CROPS, help='Crop')
    args = parser.parse_args()

    print(f"Generating global crop look-up tables for {args.crop}...")
    main(args.crop)
