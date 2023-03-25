#!usr/bin/env python3

import geopandas as gpd
import numpy as np
import pandas as pd
import rioxarray
from setting import GADM_GPKG
from setting import GADM_CSV
from setting import EARTHSTAT_CROP
from setting import FILTERED
from setting import GEOTIFFS
from setting import RAINFED_CROPLAND
from setting import IRRIGATED_CROPLAND
from setting import MOSAIC_CROPLAND
from shapely.geometry import Point
from shapely.geometry import Polygon


def gid(gid_0, gid_1, gid_2):
    '''Get GID of an administrative region

    For Level-3 region, use GID_2; for Level-2 region, use GID_1, and for Level-1 region, use GID_0
    '''
    if (gid_1 == "") and (gid_2 == ""):
        return gid_0
    elif gid_2 == "":
        return gid_1
    else:
        return gid_2


def get_grid(coord):
    GRID_SIZE = 1.0 / 24.0

    x, y = coord.xy[0][0], coord.xy[1][0]

    x0 = max(-180, x - GRID_SIZE)
    x1 = min(180, x + GRID_SIZE)
    y0 = max(-90, y - GRID_SIZE)
    y1 = min(90, y + GRID_SIZE)

    points = [(x0, y0), (x0, y1), (x1, y1), (x1, y0), (x0, y0)]

    return Polygon([[p[0], p[1]] for p in points])


def cropland_exists(rds, geometry):

    clipped = rds.rio.clip([geometry], from_disk=True)
    df = clipped[0].to_pandas()

    return np.isin(RAINFED_CROPLAND + IRRIGATED_CROPLAND + MOSAIC_CROPLAND, df.values).any()


def main():
    """Read global administrative region boundaries

    Although a Level-3 shapefile can be created, special characters cannot be rendered correctly from the shapefile.
    Therefore, reading from the gpkg file is still preferred.
    """
    print("Read gadm file")
    t = gpd.read_file(GADM_GPKG)
    t["GID"] = t.apply(lambda x: gid(x["GID_0"], x["GID_1"], x["GID_2"]), axis=1)

    print(" Create level-3 data")
    t = t.dissolve(by="GID")
    t["GID"] = t.index
    t.sort_values("UID", inplace=True)

    # Read GeoTiff file to get the grids
    print("Read GeoTiff files")
    rds = rioxarray.open_rasterio(EARTHSTAT_CROP("maize", "HarvestedAreaFraction"), masked=True)
    df = pd.DataFrame(rds[0].to_series().rename("HarvestedAreaFraction"))

    print("Convert to points")
    df["coord"] = [Point(c[1], c[0]) for c in list(df.index)]
    points = gpd.GeoDataFrame(df, geometry="coord", crs=t.crs)

    print("Find points inside administrative regions")
    points_in_polys = gpd.tools.sjoin(points, t, predicate="within", how="left")
    points_in_polys["GID"] = points_in_polys["GID"].fillna("N/A")

    print("Get grid areas")
    points_in_polys["grid"] = points_in_polys.apply(lambda x: get_grid(x["coord"]), axis=1)
    points_in_polys = gpd.GeoDataFrame(points_in_polys, geometry="grid", crs=t.crs)
    points_in_polys["GridAreaHectares"] = points_in_polys.to_crs("+proj=cea +units=m").area / 1.0E4

    print("Write to output")
    df_out = points_in_polys[[
        "GID",
        "GridAreaHectares"
    ]]
    df_out.to_csv(FILTERED, index=False)

    print("Filter out regions without cropland")
    rds = rioxarray.open_rasterio(GEOTIFFS["lc"], masked=True)
    t["cropland"] = t.apply(lambda x: cropland_exists(rds, x["geometry"]), axis=1)
    t = t[t["cropland"] == True]

    print("Calculate centroids")
    t["centroid"] = t.centroid
    t["Lat"] = t["centroid"].y
    t["Lon"] = t["centroid"].x

    print("Calculate areas")
    t["AreaKm2"] = t.to_crs("+proj=cea +units=m").area / 1.0E6

    print("Write to output")
    df_out = t[[
        "GID",
        "NAME_0",
        "NAME_1",
        "NAME_2",
        "Lat",
        "Lon",
        "AreaKm2",
    ]]
    df_out.to_csv(
        GADM_CSV,
        float_format="%.4f",
        index=False,
    )


if __name__ == "__main__":
    main()
