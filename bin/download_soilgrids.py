#!/usr/bin/env python3

import os
import rioxarray
import subprocess
from setting import GEOTIFFS
from setting import SG_VARS
from setting import SG_DEPTHS

def main():
    """Download SoilGrids 5000-m aggregated data and convert coordinates
    """

    os.makedirs("data/SoilGrids", exist_ok=True)

    for v in SG_VARS:
        path = f"data/SoilGrids/{v}"
        os.makedirs(path, exist_ok=True)

        for d in SG_DEPTHS:
            url = f"https://files.isric.org/soilgrids/latest/data_aggregated/5000m/{v}/{v}_{d}_mean_5000.tif"

            cmd = [
                "wget",
                "-c",
                url,
                "-P",
                path,
            ]

            subprocess.run(cmd)

            xds = rioxarray.open_rasterio(GEOTIFFS["SoilGrids_org"](v, d), masked=True)
            xds = xds.rio.reproject("EPSG:4326")
            xds.rio.to_raster(GEOTIFFS["SoilGrids"](v, d))


if __name__ == "__main__":
    main()
