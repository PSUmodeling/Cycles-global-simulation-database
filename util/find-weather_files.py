#!/usr/bin/env python3

import argparse
import pandas as pd
import sys
sys.path.insert(1, './bin/')
from my_funcs import read_grids
from my_funcs import find_grid


def main(params):
    ldas = params['ldas']
    lat = params['lat']
    lon = params['lon']
    fn = params['file']

    # Read meteorological reanalysis grids
    coord, mask_array = read_grids(ldas)

    if fn is not None:
        df = pd.read_csv(
            fn,
            comment='#',
            sep='\s+',
        )

        df['WEATHER'] = df.apply(
            lambda x: find_grid(ldas, float(x["LATITUDE"]), float(x["LONGITUDE"]), coord, mask_array),
            axis=1
        )

        df[['LATITUDE', 'LONGITUDE', 'WEATHER']].to_csv(
            'results.csv',
            index=False,
        )
        print(f'Results written to results.csv.')
    else:
        print(f'Weather file for {lat}x{lon} is {find_grid(ldas, lat, lon, coord, mask_array)}')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Find weather files for given locations",
    )
    parser.add_argument(
        '--ldas',
        required=True,
        type=str,
        choices=['GLDAS', 'NLDAS', 'gridMET'],
        help='Reanalysis data',
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        '--lat',
        type=float,
        help='Latitude',
    )
    parser.add_argument(
        '--lon',
        type=float,
        help='Longitude',
    )
    group.add_argument(
        '--file',
        type=str,
        help='Name of location file',
    )
    args = parser.parse_args()

    main(vars(args))
