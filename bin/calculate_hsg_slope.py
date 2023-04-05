#!/usr/bin/env python3

import geopandas as gpd
import numpy as np
import rioxarray
from rasterio.enums import Resampling
from setting import GADM_SHP
from setting import GEOTIFFS
from setting import HSG_SLOPE_CSV
from setting import RAINFED_CROPLAND
from setting import IRRIGATED_CROPLAND
from setting import MOSAIC_CROPLAND


def calculate_hsg_slope(xds_lc, xds_hsg, xds_slope, gid, geometry):
    """Get hydrologic soil group and slope for a region
    """

    print(gid)

    # If land cover data are missing for the region, return bad values
    try:
        clipped_lc = xds_lc.rio.clip([geometry], from_disk=True)
        df_lc = clipped_lc[0].to_pandas()
    except:
        return -999, -999, -999, -999

    xds = {
        "HSG": xds_hsg,
        "slope": xds_slope,
    }

    funcs = {
        "HSG": lambda x: np.bincount(x).argmax(),
        "slope": lambda x: np.median(x),
    }

    types = {
        "HSG": int,
        "slope": float,
    }

    vars = ["HSG", "slope"]

    rainfed = {}
    irrigated = {}

    for v in vars:
        # Clip geotiff using region boundary
        try:
            clipped = xds[v].rio.clip([geometry], from_disk=True)
        except:
            rainfed[v] = irrigated[v] = -999
            continue

        # Convert to dataframe
        df = clipped.rio.reproject_match(clipped_lc, resampling=Resampling.nearest)[0].to_pandas()

        # Rainfed
        x = df[df_lc.isin(RAINFED_CROPLAND)].to_numpy()
        x = x[~np.isnan(x)].astype(types[v])

        if x.size == 0:
            rainfed[v] = -999
        else:
            rainfed[v] = funcs[v](x)

        # Irrigated
        x = df[df_lc.isin(IRRIGATED_CROPLAND)].to_numpy()
        x = x[~np.isnan(x)].astype(types[v])

        if x.size == 0:
            irrigated[v] = -999
        else:
            irrigated[v] = funcs[v](x)

        # Mosaic
        if (rainfed[v] == -999) and (irrigated[v] == -999):
            x = df[df_lc.isin(MOSAIC_CROPLAND)].to_numpy()
            x = x[~np.isnan(x)].astype(types[v])

            if x.size == 0:
                y = df.to_numpy()
                y = y[~np.isnan(y)].astype(types[v])
                if y.size == 0:
                    rainfed[v] = irrigated[v] = -999
                else:
                    rainfed[v] = irrigated[v] = funcs[v](y)
            else:
                rainfed[v] = irrigated[v] = funcs[v](x)
        elif rainfed[v] == -999:
            rainfed[v] = irrigated[v]
        elif irrigated[v] == -999:
            irrigated[v] = rainfed[v]

    return rainfed["HSG"], irrigated["HSG"], rainfed["slope"], irrigated["slope"]


def main():

    xds_lc = rioxarray.open_rasterio(GEOTIFFS["lc"], masked=True)
    xds_hsg = rioxarray.open_rasterio(GEOTIFFS["HSG"], masked=True)
    xds_slope = rioxarray.open_rasterio(GEOTIFFS["slope"], masked=True)

    # Read administrative region boundaries
    print("Read gadm file")
    df = gpd.read_file(GADM_SHP)

    df[[
        "RainfedHSG",
        "IrrigatedHSG",
        "RainfedSlope",
        "IrrigatedSlope",
    ]] = df.apply(
        lambda x: calculate_hsg_slope(xds_lc, xds_hsg, xds_slope, x["GID"], x["geometry"]),
        axis=1,
        result_type="expand",
    )

    print("Write to output")
    df_out = df[[
        "GID",
        "RainfedHSG",
        "IrrigatedHSG",
        "RainfedSlope",
        "IrrigatedSlope",
    ]]
    df_out = df_out.astype({"RainfedHSG": int, "IrrigatedHSG": int})
    df_out.to_csv(
        HSG_SLOPE_CSV,
        float_format="%.2f",
        index=False,
    )


if __name__ == "__main__":
    main()
