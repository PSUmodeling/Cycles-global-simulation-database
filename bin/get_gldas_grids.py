#!/usr/bin/env python3

import pandas as pd
from setting import CROPS
from setting import LU_TYPES
from setting import LOOKUP_CSV
from setting import GLDAS_GRIDS


def main():
    """Get all GLDAS grids
    Get all GLDAS grids that are needed to simulate all crops and irrigation types
    """

    df = pd.DataFrame()
    for crop in CROPS:
        for t in LU_TYPES:
            try:
                _df = pd.read_csv(LOOKUP_CSV(crop, t), dtype=str)

                df = pd.concat(
                    [df, _df],
                    axis=0,
                    ignore_index=True,
                )
            except:
                pass

    weathers = df['Weather'].unique()

    with open(GLDAS_GRIDS, 'w') as f:
        f.write('%-15s\t%s\n' % ('LATITUDE', 'LONGITUDE'))
        for w in weathers:
            coord = w[6:-8]
            [lat, lon] = coord.split('x')
            f.write('%-15s\t' %
                (lat[:-1] if lat[-1] == 'N' else '-' + lat[:-1])
            )
            f.write('%s\n' %
                (lon[:-1] if lon[-1] == 'E' else '-' + lon[:-1])
            )


if __name__ == "__main__":
    main()
