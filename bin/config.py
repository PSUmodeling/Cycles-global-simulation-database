from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

VERSION = '4.0'

PROJECT_DIR = Path(__file__).parent.parent
DATA_DIR = PROJECT_DIR / 'data'
TEMP_DIR = PROJECT_DIR / 'temp'
LUT_DIR = PROJECT_DIR / 'crop_lut'
LUT_CSV = lambda crop, management, region: LUT_DIR / f'{crop}_{management}_{region}_lut_{VERSION}.csv'

@dataclass
class GeoSpatialDataInfoMixin:
    file_path: Callable | str
    grid_file: Path | None
    dimensions: tuple

    @property
    def di(self) -> float:
        return 360.0 / self.dimensions[1]

    @property
    def dj(self) -> float:
        return 180.0 / self.dimensions[0]

CROPGRIDS = GeoSpatialDataInfoMixin(
    file_path=lambda crop: DATA_DIR / 'CROPGRIDS' / f'CROPGRIDSv1.08_{crop}.nc',
    grid_file=DATA_DIR / 'CROPGRIDS' / f'Countries_2018.nc',
    dimensions=(3600, 7200),
)

GAEZ = GeoSpatialDataInfoMixin(
    file_path=lambda crop, management, var: DATA_DIR / f'GAEZ+_2015/GAEZAct2015_{var}_{crop}_{management}.tif',
    grid_file=None,
    dimensions=(2160, 4320),
)

ESACCI_LC = GeoSpatialDataInfoMixin(
    file_path=DATA_DIR / 'ESACCI-LC' / 'ESACCI-LC-L4-LCCS-Map-300m-P1Y-2020-v2.1.1.tif',
    grid_file=None,
    dimensions=(64800, 129600),
)

HYSOG = GeoSpatialDataInfoMixin(
    file_path=DATA_DIR / 'HYSOGs250m' / 'HYSOGs250m.tif',
    grid_file=None,
    dimensions=(67200, 172800),
)

GMTED = GeoSpatialDataInfoMixin(
    file_path=DATA_DIR / 'GMTED2010' / 'slope_1KMmd_GMTEDmd.tif',
    grid_file=None,
    dimensions=(16800, 43200),
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

# ESA CCI Land Cover classes for cropland
RAINFED_CROPLAND: frozenset[int] = frozenset({10, 11, 12})
IRRIGATED_CROPLAND: frozenset[int] = frozenset({20})
MOSAIC_CROPLAND: frozenset[int] = frozenset({30, 40})
CROPLAND_TYPES: frozenset[int] = RAINFED_CROPLAND | IRRIGATED_CROPLAND | MOSAIC_CROPLAND
