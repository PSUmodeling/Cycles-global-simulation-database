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
from my_funcs import read_grids
from my_funcs import find_grid


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
