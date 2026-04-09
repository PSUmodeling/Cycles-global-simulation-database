import cycles.weather as weather
import pandas as pd
from config import CROPS, MANAGEMENTS
from config import LUT_CSV

CONUS_WEATHER_SOURCES = ['NLDAS', 'gridMET']

def main(params):
    for crop in CROPS:
        for management in MANAGEMENTS:
            if not LUT_CSV(crop, management, 'global').exists():
                continue

            df = pd.read_csv(LUT_CSV(crop, management, 'global'))
            df = df[(df['NAME_0'] == 'United States') & (df['NAME_1'] != 'Alaska') & (df['NAME_1'] != 'Hawaii')].drop(columns=['weather'])
            if df.empty:
                continue

            locations = list(zip(df['reference_latitude'], df['reference_longitude']))
            for ldas in CONUS_WEATHER_SOURCES:
                df[f'{ldas}_weather'] = weather.find_grids(ldas, locations=locations, screen_output=False, remove_duplicates=False)

            LUT_CSV(crop, management, 'conus').parent.mkdir(parents=True, exist_ok=True)
            df[[
                'GID',
                'NAME_0',
                'NAME_1',
                'NAME_2',
                'region_area_km2',
                'crop_area_ha',
                'harvested_area_ha',
                *[f'{r}_weather' for r in CONUS_WEATHER_SOURCES],
                'soil',
                'reference_latitude',
                'reference_longitude',
            ]].to_csv(LUT_CSV(crop, management, 'conus'), float_format='%.3f')


if __name__ == '__main__':
    main()
