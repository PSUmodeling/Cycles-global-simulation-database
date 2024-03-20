import numpy as np
from netCDF4 import Dataset
from setting import MASK_FILES
from setting import NC_FIELDS
from setting import IND_I, IND_J


def read_grids(ldas):
    """Read in meteorological reanalysis grid information

    Use mask netCDF file to read in the grids, and create a land mask to filter out open water grids.
    open water grids.
    """
    # Read in grids and elevations
    with Dataset(MASK_FILES[ldas]) as nc:
        mask_array =  nc[NC_FIELDS['MASK'][ldas]][:, :] if ldas == 'gridMET' else nc[NC_FIELDS['MASK'][ldas]][0]
        mask_array = np.ma.filled(mask_array.astype(float), np.nan)
        lats, lons = np.meshgrid(nc['lat'][:], nc['lon'][:], indexing='ij')

        if ldas != 'gridMET': mask_array[mask_array == 0] = np.nan

    # Mask sea grid lat/lon as nan
    lats[np.isnan(mask_array)] = np.nan
    lons[np.isnan(mask_array)] = np.nan

    return [lats, lons], mask_array


def closest_grid(lat, lon, coord):
    """Find closest grid to an input site
    """
    lats = coord[0]
    lons = coord[1]
    dist = np.sqrt((lons - lon)**2 + (lats - lat)**2)
    closest = np.unravel_index(np.argmin(dist, axis=None), dist.shape)

    return closest


def find_grid(ldas, lat, lon, coord, mask_array):
    """Find closest land grid to an input site

    This function finds the closest unmasked grid and the closest masked grid to the specified site. By comparing the
    two grids, it will determine if the specified grid is a land point.
    """

    closest = (IND_J(ldas, lat), IND_I(ldas, lon))

    if np.isnan(mask_array[closest]):
        closest = closest_grid(lat, lon, coord)
        print(f'Nearest {ldas} land grid to {lat:.3f}x{lon:.3f} is {coord[0][closest]}x{coord[1][closest]}')

    _lat = coord[0][closest]
    _lon = coord[1][closest]

    grid = f'{ldas}_%.3f%sx%.3f%s.weather' % (
        abs(_lat),
        'S' if _lat < 0.0 else 'N',
        abs(_lon),
        'W' if _lon < 0.0 else 'E'
    )

    return grid
