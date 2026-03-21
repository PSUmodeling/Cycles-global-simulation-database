from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

VERSION = '4.0'

DATA_DIR = Path('/storage/work/yzs123/analyses/Cycles-global-simulation-database/data')
TEMP_DIR = Path('/storage/work/yzs123/analyses/Cycles-global-simulation-database/temp')
LUT_DIR = Path('/storage/work/yzs123/analyses/Cycles-global-simulation-database/crop_lut')
LUT_CSV = lambda crop, management, region: LUT_DIR / f'{crop}_{management}_{region}_lut_{VERSION}.csv'

@dataclass
class CropGridInfo:
    file_path: Callable
    grid_file: Path | None
    dimensions: tuple

    @property
    def di(self) -> float:
        return 360.0 / self.dimensions[1]

    @property
    def dj(self) -> float:
        return 180.0 / self.dimensions[0]

CROPGRIDS = CropGridInfo(
    file_path=lambda crop: DATA_DIR / 'CROPGRIDS' / f'CROPGRIDSv1.08_{crop}.nc',
    grid_file=DATA_DIR / 'CROPGRIDS' / f'Countries_2018.nc',
    dimensions=(3600, 7200),
)

GAEZ = CropGridInfo(
    file_path=lambda crop, management, var: DATA_DIR / f'GAEZ._2015/GAEZAct2015_{var}_{crop}_{management}.tif',
    grid_file=None,
    dimensions=(2160, 4320),
)

GAEZ_NAMES = {
    'bean': 'Pulses',
    'cassava': 'Cassava',
    'lentil': 'Pulses',
    'maize': 'Maize',
    'millet': 'Millet',
    'potato': 'PotatoAndSweetpotato',
    'rice': 'Rice',
    'sorghum': 'Sorghum',
    'soybean': 'Soybean',
    'sweetpotato': 'PotatoAndSweetpotato',
    'wheat': 'Wheat',
}

WGS84 = 'epsg:4326'
CEA = '+proj=cea +units=m'

