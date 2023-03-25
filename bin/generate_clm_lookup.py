#!usr/bin/env python3

import argparse
import numpy as np
import os
import pandas as pd
from netCDF4 import Dataset
from setting import VERSION
from setting import CROPS
from setting import LU_TYPES
from setting import LOOKUP_DIR
from setting import LOOKUP_CSV

LOOKUP_CLM_CSV = lambda crop, lu: f"{LOOKUP_DIR}/{crop}_{lu.lower()}_eow_lookup_{VERSION}.csv"
MASK = "data/consistent-topo-fv1.9x2.5_c130424.nc"

CLM_LA1 = -90
CLM_LO1 = 0
CLM_DI = 2.5
CLM_DJ = 180.0 / 95.0

CLM_J = lambda lat: int(round((lat - CLM_LA1) / CLM_DJ))
CLM_I = lambda lon: int(round((lon - CLM_LO1) / CLM_DI))


def read_clm_grids():
    '''Read in CLM grid information from mask file

    Use mask netCDF file to read in the grids, and create a land mask to filter out open water grids.
    '''

    # Read in grids and elevations
    with Dataset(MASK) as nc:
        mask_array = nc["LANDFRAC"][:][:]
        lats, lons = np.meshgrid(nc["lat"][:], nc["lon"][:], indexing="ij")
        mask_array[mask_array == 0] = np.nan

    # Mask sea grid lat/lon as nan
    lats[np.isnan(mask_array)] = np.nan
    lons[np.isnan(mask_array)] = np.nan

    return [lats, lons], mask_array


def closest_grid(lat, lon, coord):
    '''Find closest grid to an input site
    '''
    lats = coord[0]
    lons = coord[1]
    dist = np.sqrt((lons - lon)**2 + (lats - lat)**2)
    closest = np.unravel_index(np.argmin(dist, axis=None), dist.shape)

    return closest


def find_grid(lat, lon, coord, mask_array):
    '''Find closest land grid to an input site

    This function finds the closest unmasked grid and the closest masked grid to the specified site. By comparing the
    two grids, it will determine if the specified grid is a land point.
    '''

    closest = (CLM_J(lat), CLM_I(lon))

    if np.isnan(mask_array[closest]):   # If closest grid is a sea/non-CONUS grid
        closest = closest_grid(lat, lon, coord)
        print(f"Nearest NLDAS land grid to {lat:.3f}x{lon:.3f} is {coord[0][closest]}x{coord[1][closest]}")

    _lat = coord[0][closest]
    _lon = coord[1][closest]

    grid = f'%.1f%sx%.1f%s' % (
        abs(_lat),
        'S' if _lat < 0.0 else 'N',
        360.0 - _lon if _lon > 180.0 else _lon,
        'W' if _lon > 180.0 else 'E')

    return grid


def write_output(df, crop, lu):
    df_out = df[[
        "GID",
        "NAME_0",
        "NAME_1",
        "NAME_2",
        "Lat",
        "Lon",
        "AreaKm2",
        "AreaFraction",
        "CropLat",
        "CropLon",
        "Weather",
        "Soil",
    ]].copy()

    df_out.to_csv(
        LOOKUP_CLM_CSV(crop, lu),
        index=False,
    )


def main(crop):

    os.makedirs(LOOKUP_DIR, exist_ok=True)

    # Read NLDAS grids
    coord, mask_array = read_clm_grids()

    for t in LU_TYPES:
        # Read global look-up table
        try:
            global_df = pd.read_csv(LOOKUP_CSV(crop, t), dtype=str)
            global_df = global_df.drop(columns=["Weather"])

            global_df["Weather"] = global_df.apply(lambda x: find_grid(float(x["CropLat"]), float(x["CropLon"]), coord, mask_array), axis=1)

            write_output(global_df, crop, t)
        except:
            pass



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create crop lookup files for global administrative regions with CLM grids")
    parser.add_argument(
        "--crop",
        default="maize",
        choices=list(CROPS.keys()),
        help="Crop",
    )
    args = parser.parse_args()

    main(args.crop)
