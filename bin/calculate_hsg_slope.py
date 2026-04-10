from __future__ import annotations
import numpy as np
import pandas as pd
import rioxarray
from cycles.gadm import read_gadm
from dataclasses import dataclass
from rasterio.enums import Resampling
from shapely.geometry import Polygon
from typing import Callable
from config import ESACCI_LC, GADM, HYSOG, GMTED
from config import RAINFED_CROPLAND, IRRIGATED_CROPLAND, MOSAIC_CROPLAND
from config import HSG_SLOPE_CSV

HSG_DICT: dict[int, str] = {
    1:  'A',
    2:  'B',
    3:  'C',
    4:  'D',
    11: 'A/D',
    12: 'B/D',
    13: 'C/D',
    14: 'B/D',
}
LAND_USES: dict[str, set] = {
    'rainfed':  RAINFED_CROPLAND,
    'irrigated': IRRIGATED_CROPLAND,
    'mosaic':   MOSAIC_CROPLAND,
}
VARIABLES: tuple[str, str] = ('HSG', 'slope')
OUTPUT_COLS = ['GID', *[f'{lu}_{v}' for lu in LAND_USES for v in VARIABLES]]

@dataclass
class GeoSpatialData:
    xds: rioxarray.DataArray
    func: Callable[[np.ndarray], object] | None
    data_type: type

def parse_values(variable: str, gs_data: GeoSpatialData, reference_xds: rioxarray.DataArray, geometry: Polygon) -> dict[str: float]:
    """ Clip gs_data to geometry, reproject to match the land-cover reference, and reduce with gs_data.func.
    Returns np.nan if no data falls within the geometry or land-use mask.
    """
    try:
        clipped = gs_data.xds.rio.clip([geometry], from_disk=True)
    except Exception:
        return {f'{lu}_{variable}': np.nan for lu in LAND_USES}

    reprojected = clipped.rio.reproject_match(reference_xds, resampling=Resampling.nearest)

    results = {}
    for lu, land_use_classes in LAND_USES.items():
        lu_mask = reference_xds[0].to_pandas().isin(land_use_classes)
        values = reprojected[0].to_pandas()[lu_mask].to_numpy()
        values = values[~np.isnan(values)].astype(gs_data.data_type)
        results[f'{lu}_{variable}'] = gs_data.func(values) if values.size > 0 else np.nan

    return results


def calculate_hsg_slope(gs_datasets: dict[str, GeoSpatialData], geometry: Polygon) -> dict[str, object]:
    """ Compute hydrologic soil group and slope for each land-use type within the given geometry.
    Returns a flat dict of the form:
        {<land_use>_<variable>: value, ...}
    """
    hsg_slope = {}
    clipped_lc = gs_datasets['land_cover'].xds.rio.clip([geometry], from_disk=True)
    for variable in VARIABLES:
        hsg_slope.update(parse_values(variable, gs_datasets[variable], reference_xds=clipped_lc, geometry=geometry))

    return hsg_slope


def load_gs_datasets() -> dict[str, GeoSpatialData]:
    return {
        'land_cover': GeoSpatialData(
            xds=rioxarray.open_rasterio(ESACCI_LC.path, masked=True),
            func=None,
            data_type=int,
        ),
        'HSG': GeoSpatialData(
            xds=rioxarray.open_rasterio(HYSOG.path, masked=True),
            func=lambda x: HSG_DICT[np.bincount(x).argmax()],
            data_type=int,
        ),
        'slope': GeoSpatialData(
            xds=rioxarray.open_rasterio(GMTED.path, masked=True),
            func=np.median,
            data_type=float,
        ),
    }


def main():
    gs_datasets = load_gs_datasets()
    gadm = read_gadm(GADM.path, 'global', 'county').reset_index()

    HSG_SLOPE_CSV.parent.mkdir(parents=True, exist_ok=True)

    first_write = True
    for country, group in gadm.groupby('GID_0'):
        rows = []
        for _, row in group.iterrows():
            result = calculate_hsg_slope(gs_datasets, row['geometry'])
            rows.append({'GID': row['GID'], **result})

        country_df = pd.DataFrame(rows)[OUTPUT_COLS]
        country_df.to_csv(
            HSG_SLOPE_CSV,
            mode='w' if first_write else 'a',
            header=first_write,
            index=False,
            float_format='%.2f',
        )
        first_write = False

if __name__ == '__main__':
    main()
