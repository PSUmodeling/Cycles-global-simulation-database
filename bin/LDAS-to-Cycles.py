#!/usr/bin/env python3

import argparse
import glob
import numpy as np
import os
import subprocess
from datetime import timedelta, datetime
from netCDF4 import Dataset
from setting import DATA_DIR
from setting import GLDAS_GRIDS

WEATHER_DIR = "./weather"
COOKIE_FILE = "./.urs_cookies"

START_DATES = {
    "GLDAS": datetime.strptime("2000-01-01", "%Y-%m-%d"),
    "NLDAS": datetime.strptime("1979-01-01", "%Y-%m-%d"),
}
START_HOURS = {
    "GLDAS": 3,
    "NLDAS": 13,
}
ELEV_URLS = {
    "GLDAS": "https://ldas.gsfc.nasa.gov/sites/default/files/ldas/gldas/ELEV/GLDASp5_elevation_025d.nc4",
    "NLDAS": "https://ldas.gsfc.nasa.gov/sites/default/files/ldas/nldas/NLDAS_elevation.nc4"
}
ELEV_FILES = {
    "GLDAS": os.path.basename(ELEV_URLS['GLDAS']),
    "NLDAS": os.path.basename(ELEV_URLS['NLDAS']),
}
MASK_URL = "https://ldas.gsfc.nasa.gov/sites/default/files/ldas/nldas/NLDAS_masks-veg-soil.nc4"
MASK_FILE = "NLDAS_masks-veg-soil.nc4"
URLS = {
    "GLDAS": "https://hydro1.gesdisc.eosdis.nasa.gov/data/GLDAS/GLDAS_NOAH025_3H.2.1",
    "NLDAS": "https://hydro1.gesdisc.eosdis.nasa.gov/data/NLDAS/NLDAS_FORA0125_H.2.0",
}
EXTENSIONS = {
    "GLDAS": "nc4",
    "NLDAS": "nc",
}
NC_PREFIXES = {
    "GLDAS": "GLDAS_NOAH025_3H.A",
    "NLDAS": "NLDAS_FORA0125_H.A",
}
NC_SUFFIXES = {
    "GLDAS": "021.nc4",
    "NLDAS": "020.nc",
}
INTERVALS = {       # Data interval in hours
    "GLDAS": 3,
    "NLDAS": 1,
}
NC_FIELDS = {
    "ELEV": {
        "GLDAS": "GLDAS_elevation",
        "NLDAS": "NLDAS_elev",
    },
    "MASK": {
        "GLDAS": "",
        "NLDAS": "CONUS_mask",
    },
    "PRCP": {
        "GLDAS": "Rainf_f_tavg",
        "NLDAS": "Rainf",
    },
    "TMP": {
        "GLDAS": "Tair_f_inst",
        "NLDAS": "Tair",
    },
    "Q": {
        "GLDAS": "Qair_f_inst",
        "NLDAS": "Qair",
    },
    "UWIND": {
        "GLDAS": "Wind_f_inst",
        "NLDAS": "Wind_E",
    },
    "VWIND": {
        "GLDAS": "Wind_f_inst",
        "NLDAS": "Wind_N",
    },
    "SOLAR": {
        "GLDAS": "SWdown_f_tavg",
        "NLDAS": "SWdown",
    },
    "LONGWAVE": {
        "GLDAS": "LWdown_f_tavg",
        "NLDAS": "LWdown",
    },
    "PRES": {
        "GLDAS": "Psurf_f_inst",
        "NLDAS": "PSurf",
    },
}


def ldas_download(ldas, d):
    '''Download LDAS forcing files

    Download LDAS netCDF forcing files from GES DISC.
    '''
    print(f"    Downloading {d.strftime('%Y-%m-%d')} data...")

    # Check number of files already downloaded
    nof = len(glob.glob1(f"{DATA_DIR}/{d.strftime('%Y/%j')}", f"*.{EXTENSIONS[ldas]}"))
    # Number of files available from GES DISC
    if d == START_DATES[ldas]:
        nof_avail = int((24 - START_HOURS[ldas]) / INTERVALS[ldas])
    else:
        nof_avail = int(24 / INTERVALS[ldas])

    # If all available files are downloaded, skip
    if nof != nof_avail:
        cmd = [
            "wget",
            "--load-cookies",
            COOKIE_FILE,
            "--save-cookies",
            COOKIE_FILE,
            "--keep-session-cookies",
            "--no-check-certificate",
            "-r",
            "-c",
            "-nH",
            "-nd",
            "-np",
            "-A",
            EXTENSIONS[ldas],
            f"{URLS[ldas]}/{d.strftime('%Y/%j')}/",
            "-P",
            f"{DATA_DIR}/{d.strftime('%Y/%j')}",
        ]
        subprocess.run(
            cmd,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )


def read_nldas_grids():
    '''Read in NLDAS grid information from elevation files

    Use elevation/grid netCDF files to read in the grids and elevations. Then create a land mask to filter out open
    water grids.
    '''
    # Download mask file
    cmd = [
        "wget",
        "-N",       # Avoid downloading new copies if file already exists
        MASK_URL,
        "-P",
        DATA_DIR,
    ]
    subprocess.run(
        cmd,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

    # Download elevation file
    cmd = [
        "wget",
        "-N",       # Avoid downloading new copies if file already exists
        ELEV_URLS["NLDAS"],
        "-P",
        DATA_DIR,
    ]
    subprocess.run(
        cmd,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

    # Read in grids and elevations
    with Dataset(f"{DATA_DIR}/{ELEV_FILES['NLDAS']}") as nc:
        elev_array = nc[NC_FIELDS["ELEV"]["NLDAS"]][0]
        #elev_array = np.ma.filled(elev_array.astype(float), np.nan)

        elevs = elev_array.flatten()

    with Dataset(f"{DATA_DIR}/{MASK_FILE}") as nc:
        mask_array = nc[NC_FIELDS["MASK"]["NLDAS"]][0]

        lats, lons = np.meshgrid(nc["lat"][:], nc["lon"][:], indexing="ij")

        lats = lats.flatten()
        lons = lons.flatten()
        masks = mask_array.flatten()

    grids = []

    for k in range(mask_array.size):
        #if not math.isnan(elev_array[np.unravel_index(k, elev_array.shape)]):
        if masks[k] == 1:
            grids.append(k)

    return [lats, lons], grids, elevs


def closest_grid(lat, lon, coord):
    '''Find closest grid to an input site
    '''
    lats = coord[0]
    lons = coord[1]
    dist = np.sqrt((lons - lon)**2 + (lats - lat)**2)
    closest = np.unravel_index(np.argmin(dist, axis=None), dist.shape)[0]

    return closest


def read_gldas_grids():
    '''Read in GLDAS grids from a location file
    '''
    # Download elevation file
    cmd = [
        "wget",
        "-N",       # Avoid downloading new copies if file already exists
        ELEV_URLS["GLDAS"],
        "-P",
        DATA_DIR,
    ]
    subprocess.run(
        cmd,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

    # Read in grids and elevations
    with Dataset(f"{DATA_DIR}/{ELEV_FILES['GLDAS']}") as nc:
        elev_array = nc[NC_FIELDS["ELEV"]["GLDAS"]][0]
        elevs = elev_array.flatten()

        lats, lons = np.meshgrid(nc["lat"][:], nc["lon"][:], indexing="ij")

        lats = lats.flatten()
        lons = lons.flatten()

    grids = []

    with open(GLDAS_GRIDS) as f:
        for line in f:
            line = line.strip()

            # Skip header line, comment lines and empty lines
            if (line.startswith("#") or line.startswith("L") or (not line)):
                continue

            # Read lat/lon from location file
            strs = line.split()
            lat = float(strs[0])
            lon = float(strs[1])

            # Find the closest LDAS grid
            grid_ind = closest_grid(lat, lon, [lats, lons])

            grids.append(grid_ind)

    return [lats, lons], grids, elevs


def get_file_names(ldas, coord, grids):
    lats = coord[0]
    lons = coord[1]

    fns = []
    strs = []
    for k in grids:
        # Get lat/lon and elevation
        grid_lat = lats[k]
        grid_lon = lons[k]

        _grid = "%.3f%sx%.3f%s" % (
            abs(grid_lat), "S" if grid_lat < 0.0 else "N", abs(grid_lon), "W" if grid_lon < 0.0 else "E"
        )

        fns.append(f"{WEATHER_DIR}/{ldas}_{_grid}.weather")
        strs.append([])

    return fns, strs


def init_weather_files(coord, strs, grids, elevs):
    '''Create meteorological files and write headers
    '''
    print("Writing headers")

    # Create weather directory for meteorological files
    os.makedirs(WEATHER_DIR, exist_ok=True)

    lats = coord[0]

    # Generate meteorological files
    for kgrid in range(len(grids)):
        ## Get lat/lon and elevation
        grid_lat = lats[grids[kgrid]]
        elevation = elevs[grids[kgrid]]

        # Open meteorological file and write header lines
        strs[kgrid].append("%-23s\t%.2f\n" % ("LATITUDE", grid_lat))
        strs[kgrid].append("%-23s\t%.2f\n" % ("ALTITUDE", elevation))
        strs[kgrid].append("%-23s\t%.1f\n" % ("SCREENING_HEIGHT", 2.0))
        strs[kgrid].append("%-7s\t%-7s\t%-7s\t%-7s\t%-7s\t%-7s\t%-7s\t%-7s\t%s\n" %
            ("YEAR", "DOY", "PP", "TX", "TN", "SOLAR", "RHX", "RHN", "WIND"))
        strs[kgrid].append("%-7s\t%-7s\t%-7s\t%-7s\t%-7s\t%-7s\t%-7s\t%-7s\t%s\n" %
            ("####", "###", "mm", "degC", "degC", "MJ/m2", "%", "%", "m/s"))


def process_day(ldas, t0, grids, strs):
    '''Process one day of LDAS data and write them to meteorological files
    '''
    # Arrays to store daily values
    var = {
        "PRCP": [],
        "TMP": [],
        "WIND": [],
        "SOLAR": [],
        "LONGWAVE": [],
        "RH": [],
        "PRES": [],
    }

    print(datetime.strftime(t0, "    %Y-%m-%d"))

    t = t0
    while t < t0 + timedelta(days=1):
        if t >= START_DATES[ldas] + timedelta(hours=START_HOURS[ldas]):
            # netCDF file name
            fn = f"{t.strftime('%Y/%j')}/{NC_PREFIXES[ldas]}{t.strftime('%Y%m%d.%H%M')}.{NC_SUFFIXES[ldas]}"

            # Read one netCDF file
            with Dataset(f"{DATA_DIR}/{fn}") as nc:
                _var = read_var(ldas, grids, nc)
                for key in var:
                    var[key].append(_var[key])

        t += timedelta(hours=INTERVALS[ldas])

    # Write to Cycles weather files
    ## Process daily values
    prcp = np.array(var["PRCP"]).mean(axis=0)
    tx = np.array(var["TMP"]).max(axis=0)
    tn = np.array(var["TMP"]).min(axis=0)
    wind = np.array(var["WIND"]).mean(axis=0)
    solar = np.array(var["SOLAR"]).mean(axis=0)
    rhx = np.array(var["RH"]).max(axis=0)
    rhn = np.array(var["RH"]).min(axis=0)

    for kgrid in range(len(grids)):
        strs[kgrid].append("%-15s\t%-7s\t%-7.2f\t%-7.2f\t%-7.3f\t%-7.2f\t%-7.2f\t%-8.2f\n" % (
            t0.strftime("%Y   \t%j"),
            ("%.5f" % (prcp[kgrid] * 86400.0))[0:6],
            tx[kgrid] - 273.15,
            tn[kgrid] - 273.15,
            solar[kgrid] * 86400.0 / 1.0E6,
            rhx[kgrid] * 100.0,
            rhn[kgrid] * 100.0,
            wind[kgrid])
        )


def read_var(ldas, grids, nc):
    '''Read meteorological variables of an array of desired grids from netCDF

    The netCDF variable arrays are flattened to make reading faster
    '''
    _prcp  = nc[NC_FIELDS["PRCP"][ldas]][0].flatten()[grids]
    if ldas == "NLDAS":     # NLDAS precipitation unit is kg m-2. Convert to kg m-2 s-1 to be consistent with GLDAS
        _prcp /= INTERVALS[ldas] * 3600.0
    _tmp  = nc[NC_FIELDS["TMP"][ldas]][0].flatten()[grids]
    _uwind  = nc[NC_FIELDS["UWIND"][ldas]][0].flatten()[grids]
    _vwind  = nc[NC_FIELDS["VWIND"][ldas]][0].flatten()[grids]
    _solar = nc[NC_FIELDS["SOLAR"][ldas]][0].flatten()[grids]
    _pres  = nc[NC_FIELDS["PRES"][ldas]][0].flatten()[grids]
    _spfh  = nc[NC_FIELDS["Q"][ldas]][0].flatten()[grids]
    _longwave  = nc[NC_FIELDS["LONGWAVE"][ldas]][0].flatten()[grids]

    # Calculate relative humidity from specific humidity
    es = 611.2 * np.exp(17.67 * (_tmp - 273.15) / (_tmp - 273.15 + 243.5))
    ws = 0.622 * es / (_pres - es)
    w = _spfh / (1.0 - _spfh)
    _rh = w / ws
    _rh = np.minimum(_rh, np.full(_rh.shape, 1.0))      # Maximum RH 1.0
    _rh = np.maximum(_rh, np.full(_rh.shape, 0.01))     # Minimum RH 0.01

    _wind = np.sqrt(_uwind ** 2 + _vwind **2) if ldas == "NLDAS" else _uwind

    _var = {
        "PRCP": np.array(_prcp),            # kg m-2 s-1
        "TMP": np.array(_tmp),              # K
        "WIND": np.array(_wind),            # m s-1
        "SOLAR": np.array(_solar),          # W m-2
        "LONGWAVE": np.array(_longwave),    # W m-2
        "RH": np.array(_rh),                # -
        "PRES": np.array(_pres),            # Pa
    }

    return _var

def write_to_files(m, fns, strs):
    for k in range(len(fns)):
        with open(fns[k], mode=m) as f:
            f.writelines(strs[k])


def main(params):
    start_date = params["start"]
    end_date = params["end"]
    ldas = params["ldas"]

    # Download elevation file to get grid info
    print("Read grids")
    if ldas == "NLDAS":
        coord, grids, elev_array = read_nldas_grids()
    else:
        coord, grids, elev_array = read_gldas_grids()

    print("Initialize files")
    fns, strs = get_file_names(ldas, coord, grids)

    print(f"Create {ldas} forcing data from {datetime.strftime(start_date, '%Y-%m-%d')} to "
        f"{datetime.strftime(end_date, '%Y-%m-%d')}:")

    # Write headers if needed
    if start_date == START_DATES[ldas]:
        init_weather_files(coord, strs, grids, elev_array)

    # Download and process data day-by-day
    ## Create a cookie file. This file will be used to persist sessions across calls to Wget or Curl
    cmd = [
        "touch",
        COOKIE_FILE,
    ]
    subprocess.run(cmd)

    cday = start_date
    while cday <= end_date:
        ## Download data
        ldas_download(ldas, cday)

        ## Process each day's data
        process_day(ldas, cday, grids, strs)

        cday += timedelta(days=1)

    # Write to weather files
    write_to_files("w" if start_date == START_DATES[ldas] else "a", fns, strs)


def _main():
    parser = argparse.ArgumentParser(description="Download and parse LDAS forcing")
    parser.add_argument(
        "--start",
        default="2000-01-01",
        type=lambda s: datetime.strptime(s, '%Y-%m-%d'),
        help="Start year and month YYYY-MM-DD",
    )
    parser.add_argument(
        "--end",
        default="2000-01-01",
        type=lambda s: datetime.strptime(s, '%Y-%m-%d'),
        help="End year and month YYYY-MM-DD",
    )
    parser.add_argument(
        "--ldas",
        default="NLDAS",
        type=str,
        help="LDAS",
    )
    args = parser.parse_args()

    main(vars(args))


if __name__ == "__main__":
    _main()
