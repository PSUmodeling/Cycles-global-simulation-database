#!/usr/bin/env python3

import argparse
import numpy as np
import os
import pandas as pd
from netCDF4 import Dataset
from setting import CROPS
from setting import LU_TYPES
from setting import LOOKUP_DIR
from setting import LOOKUP_CSV
from setting import LOOKUP_CONUS_CSV
from setting import NLDAS_MASK

CROP_LA1 = 89.0 + 23.0 / 24.0
CROP_LO1 = 179.0 + 23.0 / 24.0
CROP_DI = 1.0 / 12.0
CROP_DJ = 1.0 / 12.0
CROP_NI = 4320

CROP_I = lambda x: int(round((x + CROP_LO1) / CROP_DI))
CROP_J = lambda y: int(round((CROP_LA1 - y) / CROP_DJ))
CROP_IDX = lambda x, y: CROP_J(y) * CROP_NI + CROP_I(x)

NLDAS_LA1 = 25.0625
NLDAS_LO1 = -124.9375
NLDAS_DI = 0.125
NLDAS_DJ = 0.125

NLDAS_J = lambda lat: int(round((lat - NLDAS_LA1) / NLDAS_DJ))
NLDAS_I = lambda lon: int(round((lon - NLDAS_LO1) / NLDAS_DI))


def read_nldas_grids():
    '''Read in NLDAS grid information from mask file

    Use mask netCDF file to read in the grids, and create a land mask to filter out open water grids.
    '''

    # Read in grids and elevations
    with Dataset(NLDAS_MASK) as nc:
        mask_array = nc["CONUS_mask"][0]
        mask_array = np.ma.filled(mask_array.astype(float), np.nan)
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

    closest = (NLDAS_J(lat), NLDAS_I(lon))

    if np.isnan(mask_array[closest]):   # If closest grid is a sea/non-CONUS grid
        closest = closest_grid(lat, lon, coord)
        print(f"Nearest NLDAS land grid to {lat:.3f}x{lon:.3f} is {coord[0][closest]}x{coord[1][closest]}")

    _lat = coord[0][closest]
    _lon = coord[1][closest]

    grid = f"NLDAS_%.3f%sx%.3f%s.weather" % (
        abs(_lat),
        "S" if _lat < 0.0 else "N",
        abs(_lon),
        "W" if _lon < 0.0 else "E"
    )

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
        LOOKUP_CONUS_CSV(crop, lu),
        index=False,
    )


def main(crop):

    os.makedirs(LOOKUP_DIR, exist_ok=True)

    # Read NLDAS grids
    coord, mask_array = read_nldas_grids()

    for t in LU_TYPES:
        # Read global look-up table
        try:
            global_df = pd.read_csv(LOOKUP_CSV(crop, t), dtype=str)

            conus_df = global_df[
                (global_df["NAME_0"] == "United States") &
                (global_df["NAME_1"] != "Alaska") &
                (global_df["NAME_1"] != "Hawaii")
            ].drop(columns=["Weather"])

            conus_df["Weather"] = conus_df.apply(lambda x: find_grid(float(x["CropLat"]), float(x["CropLon"]), coord, mask_array), axis=1)

            write_output(conus_df, crop, t)
        except:
            pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create crop lookup files for CONUS administrative regions")
    parser.add_argument(
        "--crop",
        default="maize",
        choices=list(CROPS.keys()),
        help="Crop",
    )
    args = parser.parse_args()

    main(args.crop)
