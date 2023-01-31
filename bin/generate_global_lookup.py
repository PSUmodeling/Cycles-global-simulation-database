#!usr/bin/env python3

import argparse
import numpy as np
import os
import pandas as pd
import rioxarray
from netCDF4 import Dataset
from setting import CROPS
from setting import LU_TYPES
from setting import LOOKUP_DIR
from setting import LOOKUP_CSV
from setting import EARTHSTAT_CROP
from setting import CROP_VARS
from setting import GADM_CSV
from setting import GEOTIFFS
from setting import FILTERED
from setting import GLDAS_MASK

CROP_LA1 = 89.0 + 23.0 / 24.0
CROP_LO1 = 179.0 + 23.0 / 24.0
CROP_DI = 1.0 / 12.0
CROP_DJ = 1.0 / 12.0
CROP_NI = 4320

CROP_I = lambda x: int(round((x + CROP_LO1) / CROP_DI))
CROP_J = lambda y: int(round((CROP_LA1 - y) / CROP_DJ))
CROP_IDX = lambda x, y: CROP_J(y) * CROP_NI + CROP_I(x)

GLDAS_LA1 = -59.875
GLDAS_LO1 = -179.875
GLDAS_DI = 0.25
GLDAS_DJ = 0.25

GLDAS_J = lambda lat: int(round((lat - GLDAS_LA1) / GLDAS_DJ))
GLDAS_I = lambda lon: int(round((lon - GLDAS_LO1) / GLDAS_DI))

_country = ""

def read_gldas_grids():
    '''Read in GLDAS grid information from mask/elevation file

    Use mask/elevation netCDF file to read in the grids, and create a land mask to filter out open water grids.
    '''

    # Read in grids and elevations
    with Dataset(GLDAS_MASK) as nc:
        mask_array = nc["GLDAS_elevation"][0]
        mask_array = np.ma.filled(mask_array.astype(float), np.nan)
        lats, lons = np.meshgrid(nc["lat"][:], nc["lon"][:], indexing="ij")

    # Mask sea grid lat/lon as nan
    lats[np.isnan(mask_array)] = np.nan
    lons[np.isnan(mask_array)] = np.nan

    return [lats, lons], mask_array


def closest_grid(lat, lon, coord):
    '''Find closest grid to an input site
    '''
    lats = coord[0]
    lons = coord[1]
    dist = np.sqrt((lons - lon)**2 + (lats - lat)**2)
    closest = np.unravel_index(np.argmin(dist, axis=None), dist.shape)

    return closest


def find_grid(lat, lon, coord, mask_array):
    '''Find closest land grid to an input site

    This function finds the closest unmasked grid and the closest masked grid to the specified site. By comparing the
    two grids, it will determine if the specified grid is a land point.
    '''

    closest = (GLDAS_J(lat), GLDAS_I(lon))

    if np.isnan(mask_array[closest]):   # If closest grid is a sea/non-CONUS grid
        closest = closest_grid(lat, lon, coord)
        print(f"Nearest GLDAS land grid to {lat:.3f}x{lon:.3f} is {coord[0][closest]}x{coord[1][closest]}")

    _lat = coord[0][closest]
    _lon = coord[1][closest]

    grid = f"GLDAS_%.3f%sx%.3f%s.weather" % (
        abs(_lat),
        "S" if _lat < 0.0 else "N",
        abs(_lon),
        "W" if _lon < 0.0 else "E"
    )

    return grid


def calculate_crop_fractions(points, crop, coord, mask_array, country, gid, lat, lon):

    global _country

    if (country != _country):
        print(country)
        _country = country

    # Default values
    fractions = {}
    lats = {}
    lons = {}
    grids = {}
    soils = {}

    for t in LU_TYPES:
        fractions[t] = 0.0
        lats[t] = -999
        lons[t] = -999
        grids[t] = ""
        soils[t] = ""

    # Find crop grids inside region
    _points = points[points["GID"] == gid]

    no_grids = False
    # Calculate harvested area fractions
    if _points.empty:   # No crop grids inside region
        no_grids = True
        _points = points.iloc[[CROP_IDX(lon, lat)]]     # Find nearest crop grid

    for t in LU_TYPES:
        if _points[f"{t}AreaFraction"].notna().values.any():
            fractions[t] = np.nanmean(_points[f"{t}AreaFraction"])

            if fractions[t] > 0.0:
                if no_grids == True:
                    lats[t], lons[t] = lat, lon
                else:
                    lats[t], lons[t] = _points[f"_{t}Production"].index[np.nanargmax(_points[f"_{t}Production"].values)]

    # Set weather and soil files
    for t in LU_TYPES:
        if fractions[t] > 0.0:
            grids[t] = find_grid(lats[t], lons[t], coord, mask_array)
            soils[t] = f"{crop}_{t.lower()}_sg_{gid}.soil"

    return (
        fractions["Rainfed"],
        fractions["Irrigated"],
        lats["Rainfed"], lons["Rainfed"],
        lats["Irrigated"], lons["Irrigated"],
        grids["Rainfed"],
        grids["Irrigated"],
        soils["Rainfed"],
        soils["Irrigated"],
    )


def write_output(df, crop, lu):
    df_out = df[[
        "GID",
        "NAME_0",
        "NAME_1",
        "NAME_2",
        "Lat",
        "Lon",
        "AreaKm2",
        f"{lu}AreaFraction",
        f"{lu}CropLat",
        f"{lu}CropLon",
        f"{lu}Weather",
        f"{lu}Soil",
    ]].copy()

    df_out = df_out.rename(columns=lambda x: x[len(lu):] if x.startswith(lu) else x)

    df_out = df_out[df_out["AreaFraction"] > 0.0]
    if not df_out.empty:
        df_out["AreaFraction"] = df_out["AreaFraction"].map(lambda x: "%.6g" % x)

        df_out.to_csv(
            LOOKUP_CSV(crop, lu),
            float_format="%.4f",
            index=False,
        )


def main(crop):

    os.makedirs(LOOKUP_DIR, exist_ok=True)

    # Read GLDAS grids
    coord, mask_array = read_gldas_grids()

    # Read administrative region boundaries
    print("Read gadm file")
    gadm_df = pd.read_csv(GADM_CSV)

    # Read filtered crop GeoTIff grids
    filtered = pd.read_csv(FILTERED, dtype={"GID": "str"})

    # EarthStat crop harvested fraction and yield database GeoTiff images are not consistent. Some GeoTiff images label
    # water grids as NaN while others don't. To be consistent, read in the GeoTiff image of maize, which has water grids
    # corrected labeled as NaN, to mark other crop water grids NaN.
    ref_rds = rioxarray.open_rasterio(EARTHSTAT_CROP("maize", "HarvestedAreaHectares"), masked=True)
    ref_df = ref_rds[0].to_series()

    # Read all GeoTiff files and merge into one DataFrame
    print("Read GeoTiff files")
    rds = {}
    for v in CROP_VARS:
        rds[v] = rioxarray.open_rasterio(GEOTIFFS[v](crop), masked=True)

        if v == CROP_VARS[0]:
            crop_df = pd.DataFrame(rds[v][0].to_series().rename(v))
        else:
            rds[v] = rds[v].rio.reproject_match(rds[CROP_VARS[0]])  # Project GeoTiff grids to the same grids to assure
                                                                    # correct reading of data
            crop_df = crop_df.join(pd.DataFrame(rds[v][0].to_series().rename(v)))

    # Label water grids NaN
    crop_df[ref_df.isna()] = np.nan

    # Add region and grid area information
    crop_df["GID"] = filtered["GID"].to_list()
    #crop_df["GridAreaHectares"] = filtered["GridAreaHectares"].to_list()

    # Calculate rainfed and irrigated area fractions
    crop_df["_RainfedAreaOverTotal"] = crop_df["HarvestedAreaRainfed"] / (crop_df["HarvestedAreaRainfed"] + crop_df["HarvestedAreaIrrigated"])
    crop_df["_RainfedAreaOverTotal"] = crop_df["_RainfedAreaOverTotal"].fillna(1.0)

    crop_df["RainfedAreaFraction"] = crop_df["HarvestedAreaFraction"] * crop_df["_RainfedAreaOverTotal"]
    crop_df["IrrigatedAreaFraction"] = crop_df["HarvestedAreaFraction"] - crop_df["RainfedAreaFraction"]

    crop_df["_RainfedProductionFraction"] = crop_df["ProductionRainfed"] / (crop_df["ProductionRainfed"] + crop_df["ProductionIrrigated"])
    crop_df["_RainfedProductionFraction"] = crop_df["_RainfedProductionFraction"].fillna(1.0)

    crop_df["_RainfedProduction"] = crop_df["Production"] * crop_df["_RainfedProductionFraction"]
    crop_df["_IrrigatedProduction"] = crop_df["Production"] - crop_df["_RainfedProduction"]

    #crop_df = crop_df.drop(CROP_VARS, axis=1)

    print("Calculate harvested area fractions")
    gadm_df[[
        "RainfedAreaFraction",
        "IrrigatedAreaFraction",
        "RainfedCropLat",
        "RainfedCropLon",
        "IrrigatedCropLat",
        "IrrigatedCropLon",
        "RainfedWeather",
        "IrrigatedWeather",
        "RainfedSoil",
        "IrrigatedSoil",
    ]] = gadm_df.apply(
        lambda x: calculate_crop_fractions(crop_df, crop, coord, mask_array, x["NAME_0"], x["GID"], x["Lat"], x["Lon"]),
        axis=1,
        result_type="expand",
    )

    # Write to csv files
    print("Write to output")
    for t in LU_TYPES:
        write_output(gadm_df, crop, t)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create crop lookup files for global administrative regions")
    parser.add_argument(
        "--crop",
        default="maize",
        choices=list(CROPS.keys()),
        help="Crop",
    )
    args = parser.parse_args()

    main(args.crop)

