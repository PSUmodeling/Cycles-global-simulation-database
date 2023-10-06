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
from setting import IND_I, IND_J
from setting import MASK_FILES
from setting import NC_FIELDS


def read_grids(ldas):
    """Read in meteorological reanalysis grid information

    Use mask netCDF file to read in the grids, and create a land mask to filter out open water grids.
    open water grids.
    """
    # Read in grids and elevations
    with Dataset(MASK_FILES[ldas]) as nc:
        mask_array = nc[NC_FIELDS['MASK'][ldas]][0] if ldas == 'NLDAS' else nc[NC_FIELDS['MASK'][ldas]][:, :]
        mask_array = np.ma.filled(mask_array.astype(float), np.nan)
        lats, lons = np.meshgrid(nc['lat'][:], nc['lon'][:], indexing='ij')

        if ldas == 'NLDAS': mask_array[mask_array == 0] = np.nan

    # Mask sea grid lat/lon as nan
    lats[np.isnan(mask_array)] = np.nan
    lons[np.isnan(mask_array)] = np.nan

    return [lats, lons], mask_array


def closest_grid(lat, lon, coord):
    """Find closest grid to an input site
    """
    lats = coord[0]
    lons = coord[1]
    dist = np.sqrt((lons - lon)**2 + (lats - lat)**2)
    closest = np.unravel_index(np.argmin(dist, axis=None), dist.shape)

    return closest


def find_grid(ldas, lat, lon, coord, mask_array):
    """Find closest land grid to an input site

    This function finds the closest unmasked grid and the closest masked grid to the specified site. By comparing the
    two grids, it will determine if the specified grid is a land point.
    """

    closest = (IND_J(ldas, lat), IND_I(ldas, lon))

    if np.isnan(mask_array[closest]):
        closest = closest_grid(lat, lon, coord)
        print(f'Nearest {ldas} land grid to {lat:.3f}x{lon:.3f} is {coord[0][closest]}x{coord[1][closest]}')

    _lat = coord[0][closest]
    _lon = coord[1][closest]

    grid = f'{ldas}_%.3f%sx%.3f%s.weather' % (
        abs(_lat),
        'S' if _lat < 0.0 else 'N',
        abs(_lon),
        'W' if _lon < 0.0 else 'E'
    )

    return grid


def write_output(ldas, df, crop, lu):
    df_out = df[[
        'GID',
        'NAME_0',
        'NAME_1',
        'NAME_2',
        'Lat',
        'Lon',
        'AreaKm2',
        'AreaFraction',
        'CropLat',
        'CropLon',
        'Weather',
        'Soil',
    ]].copy()

    df_out.to_csv(
        LOOKUP_CONUS_CSV(crop, lu, ldas),
        index=False,
    )


def main(params):
    ldas = params['ldas']

    os.makedirs(LOOKUP_DIR + '/conus', exist_ok=True)

    # Read meteorological reanalysis grids
    coord, mask_array = read_grids(ldas)

    for crop in CROPS:
        for t in LU_TYPES:
            # Read global look-up table
            try:
                global_df = pd.read_csv(LOOKUP_CSV(crop, t), dtype=str)
                print(t, crop)

                conus_df = global_df[
                    (global_df["NAME_0"] == "United States") &
                    (global_df["NAME_1"] != "Alaska") &
                    (global_df["NAME_1"] != "Hawaii")
                ].drop(columns=["Weather"])

                conus_df["Weather"] = conus_df.apply(
                    lambda x: find_grid(ldas, float(x["CropLat"]), float(x["CropLon"]), coord, mask_array), axis=1
                )

                write_output(ldas, conus_df, crop, t)
            except:
                pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create crop lookup files for CONUS administrative regions using gridMET grids",
    )
    parser.add_argument(
        '--ldas',
        default='NLDAS',
        type=str,
        choices=['NLDAS', 'gridMET'],
        help='Reanalysis data',
    )
    args = parser.parse_args()

    main(vars(args))
