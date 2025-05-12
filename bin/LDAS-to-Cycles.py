#!/usr/bin/env python3

import argparse
import os
import cycles.weather as weather
from datetime import datetime
from setting import GRIDS
from setting import DATA_DIR

WEATHER_DIR = './weather'
COOKIE_FILE = './.urs_cookies'

def read_grids_from_file(grid_fn):
    with open(grid_fn) as f:
        strs = [line.strip() for line in f if line.strip() and line.strip()[0] != '#']

    locations = [
        (float(line.split('\t')[0]), float(line.split('\t')[1])) for line in strs[1:]
    ]

    return locations


def main(params):
    start_date = params['start']
    end_date = params['end']
    reanalysis = params['reanalysis']
    model = params['model'] if reanalysis == 'MACA' else None
    rcp = params['rcp'] if reanalysis == 'MACA' else None

    if start_date < weather.START_DATES[reanalysis]: exit('Start date error.')

    # Create weather directory for meteorological files
    os.makedirs(WEATHER_DIR, exist_ok=True)

    print('Read grids')
    if reanalysis != 'NLDAS': locations = read_grids_from_file(GRIDS(reanalysis))

    print('Download data')
    weather.download_gridmet(DATA_DIR, start_date.year, end_date.year)

    print('Process weather files')
    weather.process_gridmet(DATA_DIR, WEATHER_DIR, start_date, end_date, locations, header=False)


def _main():
    parser = argparse.ArgumentParser(description='Download and parse reanalysis forcing data')
    parser.add_argument(
        '--start',
        default='2000-01-01',
        type=lambda s: datetime.strptime(s, '%Y-%m-%d'),
        help='Start year and month YYYY-MM-DD',
    )
    parser.add_argument(
        '--end',
        default='2000-01-01',
        type=lambda s: datetime.strptime(s, '%Y-%m-%d'),
        help='End year and month YYYY-MM-DD',
    )
    parser.add_argument(
        '--reanalysis',
        default='NLDAS',
        type=str,
        choices=['GLDAS', 'NLDAS', 'gridMET', 'MACA'],
        help='Reanalysis data',
    )
    parser.add_argument(
        '--model',
        type=str,
        choices=[
            'bcc-csm1-1',
            'bcc-csm1-1-m',
            'BNU-ESM',
            'CanESM2',
            'CCSM4',
            'CNRM-CM5',
            'CSIRO-Mk3-6-0',
            'GFDL-ESM2G',
            'GFDL-ESM2M',
            'HadGEM2-CC365',
            'HadGEM2-ES365',
            'inmcm4',
            'IPSL-CM5A-LR',
            'IPSL-CM5A-MR',
            'IPSL-CM5B-LR',
            'MIROC5',
            'MIROC-ESM',
            'MIROC-ESM-CHEM',
            'MRI-CGCM3',
            'NorESM1-M',
        ],
        help='CMIP5 model',
    )
    parser.add_argument(
        '--rcp',
        type=str,
        choices=['45', '85'],
        help='RCP',
    )
    args = parser.parse_args()

    main(vars(args))


if __name__ == '__main__':
    _main()
