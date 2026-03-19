import geopandas as gpd
import numpy as np
import pandas as pd
import xarray
from pathlib import Path
from shapely.geometry import Polygon
from config import CROPGRIDS
from config import GADM_PATH
from config import WGS84, CEA

def main():
    # Load CROPGRIDS grid information and create a GeoDataFrame of grid cell polygons
    with xarray.open_dataset(CROPGRIDS.grid_file) as nc:
        lons, lats = np.meshgrid(nc['lon'].values, nc['lat'].values)

    cropgrids_df = pd.DataFrame({
        'lon': lons.flatten(),
        'lat': lats.flatten(),
        'grid_index': range(CROPGRIDS.dimensions[0] * CROPGRIDS.dimensions[1]),
        'country': nc['country'].values.flatten().astype(int)
    })
    cropgrids_df = cropgrids_df[cropgrids_df['country'] > -1].reset_index()

    cropgrids_df['geometry'] = cropgrids_df.apply(
        lambda row: Polygon([
            (row['lon'] - 0.5 * CROPGRIDS.di, row['lat'] - 0.5 * CROPGRIDS.dj),
            (row['lon'] + 0.5 * CROPGRIDS.di, row['lat'] - 0.5 * CROPGRIDS.dj),
            (row['lon'] + 0.5 * CROPGRIDS.di, row['lat'] + 0.5 * CROPGRIDS.dj),
            (row['lon'] - 0.5 * CROPGRIDS.di, row['lat'] + 0.5 * CROPGRIDS.dj)
        ]),
        axis=1,
    )

    cropgrids_gdf = gpd.GeoDataFrame(
        cropgrids_df,
        geometry=cropgrids_df['geometry'],
        crs=WGS84,
    )
    cropgrids_gdf['grid_area'] = cropgrids_gdf.to_crs(CEA).area

    gadm = gpd.read_file(GADM_PATH)
    for gid in gadm['GID_0'].unique():
        # Perform spatial intersection between CROPGRIDS grid cells and GADM regions
        gdf = gpd.overlay(gadm[gadm['GID_0'] == gid], cropgrids_gdf, how='intersection')
        gdf['grid_fraction'] = gdf.to_crs(CEA).area / gdf['grid_area']
        gdf = gdf[gdf['grid_fraction'] > 1.0E-3]

        gdf[['GID', 'grid_index', 'country', 'grid_fraction']].to_csv(
            Path('./temp') / f'{gid}.csv',
            float_format='%.3f',
            index=False,
        )


if __name__ == '__main__':
    main()
