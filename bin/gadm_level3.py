import geopandas as gpd
from config import GADM

GADM_GPKG = GADM.path / 'gadm_410.gpkg'
GADM_SHP = GADM.path / 'gadm41_global_2.shp'

def choose_gid(gid_0: str, gid_1: str, gid_2: str) -> str:
    return gid_2 or gid_1 or gid_0


def main():
    """Generate shapefile for global level-3 (or above) administrative regions
    """
    # Read global administrative region boundaries
    print("Read gadm file")
    t = gpd.read_file(GADM_GPKG)

    # Fix GID_1 and GID_2 for Ghana (GHA)
    t.loc[t['GID_0'] == 'GHA', 'GID_1'] = t[t['GID_0'] == 'GHA']['GID_1'].map(lambda x: x.replace('GHA', 'GHA.'))
    t.loc[t['GID_0'] == 'GHA', 'GID_2'] = t[t['GID_0'] == 'GHA']['GID_2'].map(lambda x: x.replace('GHA', 'GHA.'))

    # Choose GID for each Level-3 (or above) region
    t['GID'] = t.apply(lambda x: choose_gid(x['GID_0'], x['GID_1'], x['GID_2']), axis=1)

    print("Create level-3 data")
    t = t.dissolve(by='GID')
    t.sort_values('UID', inplace=True)

    t.to_file(GADM_SHP)


if __name__ == '__main__':
    main()
