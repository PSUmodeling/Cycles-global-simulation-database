from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

VERSION = '4.0'

PROJECT_DIR = Path(__file__).parent.parent
DATA_DIR = PROJECT_DIR / 'data'
TEMP_DIR = PROJECT_DIR / 'temp'
LUT_DIR = PROJECT_DIR / 'crop_lut'
LUT_CSV = lambda crop, management, region: LUT_DIR / f'{crop}_{management}_{region}_lut_{VERSION}.csv'

CROPS: frozenset(str) = frozenset({
    'bean',
    'cassava',
    'lentil',
    'maize',
    'millet',
    'potato',
    'rice',
    'sorghum',
    'soybean',
    'sweetpotato',
    'wheat',
})
MANAGEMENTS: frozenset[str] = frozenset({'rainfed', 'irrigated'})

SOIL_SOURCE = 'SoilGrids'
WEATHER_SOURCE = 'GLDAS'

@dataclass
class CropGridsData:
    path: Callable
    grid_file: Path
    dimensions: tuple=(3600, 7200)
    ocean_value: int=-1

    @property
    def di(self) -> float:
        return 360.0 / self.dimensions[1]

    @property
    def dj(self) -> float:
        return 180.0 / self.dimensions[0]

@dataclass
class GaezData:
    path: Callable
    dimensions: tuple=(2160, 4320)
    variables: dict[str, str] = {
        'harvested_area': 'HarvArea',
        'production': 'Production',
    }
    crop_names: dict[str, str] = {
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

@dataclass
class GenericRasterData:
    path: Path
    dimensions: tuple

@dataclass
class GenericData:
    path: Path

CROPGRIDS = CropGridsData(
    path=lambda crop: DATA_DIR / 'CROPGRIDS' / f'CROPGRIDSv1.08_{crop}.nc',
    grid_file=DATA_DIR / 'CROPGRIDS' / f'Countries_2018.nc',
)

GAEZ = GaezData(
    path=lambda crop, management, var: DATA_DIR / f'GAEZ+_2015/GAEZAct2015_{var}_{crop}_{management}.tif',
)

ESACCI_LC = GenericRasterData(
    path=DATA_DIR / 'ESACCI-LC' / 'ESACCI-LC-L4-LCCS-Map-300m-P1Y-2020-v2.1.1.tif',
    dimensions=(64800, 129600),
)

HYSOG = GenericRasterData(
    path=DATA_DIR / 'HYSOGs250m' / 'HYSOGs250m.tif',
    dimensions=(67200, 172800),
)

GMTED = GenericRasterData(
    path=DATA_DIR / 'GMTED2010' / 'slope_1KMmd_GMTEDmd.tif',
    dimensions=(16800, 43200),
)

GADM = GenericData(
    path=DATA_DIR / 'GADM',
)

WGS84 = 'epsg:4326'
CEA = '+proj=cea +units=m'

# ESA CCI Land Cover classes for cropland
RAINFED_CROPLAND: frozenset[int] = frozenset({10, 11, 12})
IRRIGATED_CROPLAND: frozenset[int] = frozenset({20})
MOSAIC_CROPLAND: frozenset[int] = frozenset({30, 40})
CROPLAND_TYPES: frozenset[int] = RAINFED_CROPLAND | IRRIGATED_CROPLAND | MOSAIC_CROPLAND
