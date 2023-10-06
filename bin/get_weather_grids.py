#!/usr/bin/env python3

import argparse
import pandas as pd
from setting import CROPS
from setting import LU_TYPES
from setting import LOOKUP_CSV
from setting import LOOKUP_CONUS_CSV
from setting import GRIDS


def main(params):
    """Get weather grids
    Get all GLDAS or gridMET grids that are needed to simulate all crops and irrigation types
    """

    ldas = params['ldas']

    df = pd.DataFrame()
    for crop in CROPS:
        for t in LU_TYPES:
            try:
                if ldas == 'GLDAS':
                    csv = LOOKUP_CSV(crop, t)
                elif ldas == 'gridMET':
                    csv = LOOKUP_CONUS_CSV(crop, t, ldas)

                _df = pd.read_csv(csv, dtype=str)

                df = pd.concat(
                    [df, _df],
                    axis=0,
                    ignore_index=True,
                )
            except:
                pass

    weathers = df['Weather'].unique()

    start = len(ldas) + 1

    with open(GRIDS(ldas), 'w') as f:
        f.write('%-15s\t%s\n' % ('LATITUDE', 'LONGITUDE'))
        for w in weathers:
            coord = w[start:-8]
            [lat, lon] = coord.split('x')
            f.write('%-15s\t' %
                (lat[:-1] if lat[-1] == 'N' else '-' + lat[:-1])
            )
            f.write('%s\n' %
                (lon[:-1] if lon[-1] == 'E' else '-' + lon[:-1])
            )


def _main():
    parser = argparse.ArgumentParser(description='Get weather grids')
    parser.add_argument(
        '--ldas',
        default='GLDAS',
        type=str,
        choices=['GLDAS', 'gridMET'],
        help='Reanalysis data',
    )
    args = parser.parse_args()

    main(vars(args))


if __name__ == '__main__':
    _main()
