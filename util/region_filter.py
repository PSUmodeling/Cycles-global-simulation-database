#!usr/bin/env python3

import argparse
import numpy as np
import os
import pandas as pd
import subprocess
from setting import CROPS
from setting import LU_TYPES
from setting import LOOKUP_DIR
from setting import LOOKUP_CSV
from setting import LOOKUP_CONUS_CSV
from setting import SEVEN_ZIP
from setting import SOIL_ARCHIVE

WEATHER_DIR = "weather"


def main(args):

    crop = args.crop
    range = args.range

    strs = ""
    region = []
    for s in args.region:
        strs += " " + s
    region_strs = strs.replace(",", "_").replace(" ", "")
    for s in strs.split(","):
        region.append(s.strip())

    print(f"Target region: {', '.join(region)}")

    os.makedirs(f"{region_strs}_{crop}", exist_ok=True)

    weather = []

    # Generate look-up tables
    for lu in LU_TYPES:
        csv = LOOKUP_CSV(crop, lu) if range == "global" else LOOKUP_CONUS_CSV(crop, lu)

        # Read lookup table
        try:
            df = pd.read_csv(csv)
        except:
            print(f"Error reading {csv}.")
            continue

        print(f"Filter {lu.lower()} look-up table")
        df = df[df["NAME_0"] == region[0]]

        if len(region) > 1:
            df = df[df["NAME_1"] == region[1]]
        if len(region) > 2:
            df = df[df["NAME_2"] == region[2]]

        if df.empty:
            print(f"No {lu.lower()} {crop} records found for {', '.join(region)}.")
            continue

        df.to_csv(
            f"{region_strs}_{crop}/{csv.replace(range, region_strs)[len(LOOKUP_DIR):]}",
            index=False,
        )

        # Get soil data
        print(f"Generate {lu.lower()} soil file archive")
        for w in df["Soil"]:
            cmd = [
                SEVEN_ZIP,
                "e",
                SOIL_ARCHIVE(crop, lu),
                f"-o{region_strs}_{crop}",
                w,
                "-y",
            ]
            subprocess.run(
                cmd,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )

        if len(region) != 3:
            cmd = f"{SEVEN_ZIP} a -sdel {region_strs}_{crop}/{crop}_{lu.lower()}_{region_strs}_soil.7z ./{region_strs}_{crop}/{crop}_{lu.lower()}*.soil"
            subprocess.run(
                cmd,
                shell=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )

        # Get weather data
        for w in df["Weather"].unique():
            weather.append(w)

    weather_archive = "GLDAS_2000-2021.7z" if range == "global" else "NLDAS_CONUS_1979-2022.7z"
    ldas = "GLDAS" if range == "global" else "NLDAS"

    print(f"Generate weather file archive")
    if len(weather):
        for w in np.unique(np.array(weather)):
            cmd = [
                SEVEN_ZIP,
                "e",
                f"{WEATHER_DIR}/{weather_archive}",
                f"-o{region_strs}_{crop}",
                w,
                "-y",
            ]
            subprocess.run(
                cmd,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )


        if len(region) != 3:
            cmd = f"{SEVEN_ZIP} a -sdel {region_strs}_{crop}/{ldas}_{crop}_{region_strs}_weather.7z ./{region_strs}_{crop}/*.weather"
            subprocess.run(
                cmd,
                shell=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter lookup tables using countries/regions")
    parser.add_argument(
        "--region",
        nargs="+",
        required=True,
        help="Region",
    )
    parser.add_argument(
        "--crop",
        required=True,
        choices=list(CROPS.keys()),
        help="Crop",
    )
    parser.add_argument(
        "--range",
        default="global",
        choices=["global", "conus"],
        help="Range (global or conus)",
    )
    args = parser.parse_args()

    main(args)

