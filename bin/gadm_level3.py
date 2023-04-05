#!/usr/bin/env python3

import geopandas as gpd
from setting import GADM_GPKG
from setting import GADM_SHP


def gid(gid_0, gid_1, gid_2):
    if (gid_1 == "") and (gid_2 == ""):
        return gid_0
    elif gid_2 == "":
        return gid_1
    else:
        return gid_2


def main():
    """Generate shapefile for global level-3 administrative regions
    """
    # Read global administrative region boundaries
    print("Read gadm file")
    t = gpd.read_file(GADM_GPKG)

    # Define GID for each Level-3 region
    t["GID"] = t.apply(lambda x: gid(x["GID_0"], x["GID_1"], x["GID_2"]), axis=1)

    print(" Create level-3 data")
    t = t.dissolve(by="GID")
    t.sort_values("UID", inplace=True)

    t.to_file(GADM_SHP)


if __name__ == "__main__":
    main()
