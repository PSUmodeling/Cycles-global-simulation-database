import geopandas as gpd
import numpy as np
import pandas as pd
import rioxarray
import xarray
from cycles.gadm import read_gadm
from pathlib import Path
from shapely.geometry import Polygon
from config import DATA_DIR, TEMP_DIR
from config import CROPGRIDS, ESACCI_LC, GADM
from config import CROPLAND_TYPES
from config import WGS84, CEA

GRID_FRACTION_THRESHOLD: float = 1E-6

def make_cell_polygon(lon: float, lat: float) -> Polygon:
    half_di = 0.5 * CROPGRIDS.di
    half_dj = 0.5 * CROPGRIDS.dj
    return Polygon([
        (lon - half_di, lat - half_dj),
        (lon + half_di, lat - half_dj),
        (lon + half_di, lat + half_dj),
        (lon - half_di, lat + half_dj),
    ])


def build_cropgrids_gdf() -> gpd.GeoDataFrame:
    """Load CROPGRIDS grid info and return a GeoDataFrame of cell polygons.
    """
    with xarray.open_dataset(CROPGRIDS.grid_file) as nc:
        lons, lats = np.meshgrid(nc['lon'].values, nc['lat'].values)
        country = nc['country'].values.flatten().astype(int)

    df = pd.DataFrame({
        'lon': lons.flatten(),
        'lat': lats.flatten(),
        'grid_index': range(CROPGRIDS.dimensions[0] * CROPGRIDS.dimensions[1]),
        'country': country,
    })
    df = df[df['country'] != CROPGRIDS.ocean_value].reset_index(drop=True)
    df['geometry'] = df.apply(lambda row: make_cell_polygon(row['lon'], row['lat']), axis=1)

    gdf = gpd.GeoDataFrame(df, geometry='geometry', crs=WGS84)
    gdf['grid_area'] = gdf.to_crs(CEA).area
    return gdf


def has_cropland(rds: rioxarray.DataArray, geometry: Polygon) -> bool:
    """Return True if any cropland type exists within the geometry.
    """
    clipped = rds.rio.clip([geometry], from_disk=True)
    return bool(CROPLAND_TYPES & set(clipped[0].values.flat))


def load_cropland_gadm() -> gpd.GeoDataFrame:
    """Load county-level GADM regions filtered to those containing cropland.
    """
    gadm = read_gadm(GADM.path, 'global', 'county').drop('ATA').reset_index()   # Remove Antarctica to accelerate processing
    with rioxarray.open_rasterio(ESACCI_LC.path, masked=True) as rds:
        gadm['cropland'] = gadm['geometry'].apply(lambda geom: has_cropland(rds, geom))
    return gadm[gadm['cropland']].reset_index(drop=True)


def write_country_intersections(gadm: gpd.GeoDataFrame, cropgrids_gdf: gpd.GeoDataFrame, out_dir: Path) -> None:
    """Intersect GADM regions with CROPGRIDS cells and write one CSV per country.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    for gid in gadm['GID_0'].unique():
        gdf = gpd.overlay(gadm[gadm['GID_0'] == gid], cropgrids_gdf, how='intersection')
        gdf['grid_fraction'] = gdf.to_crs(CEA).area / gdf['grid_area']
        gdf = gdf[gdf['grid_fraction'] > GRID_FRACTION_THRESHOLD]
        if gdf.empty:
            continue
        gdf[['GID', 'grid_index', 'country', 'grid_fraction']].to_csv(
            out_dir / f'{gid}.csv',
            float_format='%.6f',
            index=False,
        )


def main() -> None:
    cropgrids_gdf = build_cropgrids_gdf()
    gadm = load_cropland_gadm()
    write_country_intersections(gadm, cropgrids_gdf, TEMP_DIR)


if __name__ == '__main__':
    main()
