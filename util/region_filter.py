#!/usr/bin/env python3

import argparse
import os
import pandas as pd
import subprocess
import sys
sys.path.insert(1, './bin/')
from setting import CROPS
from setting import LU_TYPES
from setting import LOOKUP_DIR
from setting import LOOKUP_CSV
from setting import LOOKUP_CONUS_CSV
from setting import SEVEN_ZIP
from setting import SOIL_ARCHIVE
from setting import WEATHER_DIR
from setting import WEATHER_FILE_ARCHIVE


def main(args):

    crop = args.crop
    ldas = args.ldas
    directory = args.directory

    os.makedirs(directory, exist_ok=True)

    # Read regions from an input file
    lines = []
    with open('regions.csv') as f:
        _lines = f.readlines()
    lines = [line for line in _lines if not line.strip().startswith("#")]

    dfs = {}

    for lu in LU_TYPES:
        # Read lookup table
        csv = LOOKUP_CSV(crop, lu) if ldas == "GLDAS" else LOOKUP_CONUS_CSV(crop, lu, ldas)
        try:
            lookup = pd.read_csv(csv)
        except:
            print(f"Error reading {csv}.")
            continue

        print(f"Filter {lu.lower()} look-up table")

        dfs[lu] = pd.DataFrame()
        for line in lines:
            region = line.strip().split(',')

            _df = lookup[lookup["NAME_0"] == region[-1].strip()].copy()
            if len(region) > 1:
                _df = _df[_df["NAME_1"] == region[-2].strip()]
            if len(region) > 2:
                _df = _df[_df["NAME_2"] == region[-3].strip()]

            if _df.empty:
                print(f"No {lu.lower()} {crop} records found for {','.join(region)}.")
                continue

            _df.loc[:, 'IrrigationType'] = lu
            dfs[lu] = pd.concat([dfs[lu], _df])

        if dfs[lu].empty: continue

        dfs[lu] = dfs[lu].drop_duplicates()
        print(dfs[lu])

        # Get soil data
        print(f"Extract {lu.lower()} soil file archive")

        cmd = [
            SEVEN_ZIP,
            "e",
            SOIL_ARCHIVE(crop, lu),
            f"-o{directory}",
        ]

        for s in dfs[lu]["Soil"]:
            cmd.append(s)

        cmd.append("-y")

        subprocess.run(
            cmd,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )

    df = pd.DataFrame()
    for lu in LU_TYPES:
        df = pd.concat([df, dfs[lu]])

    if df.empty:
        exit('No regions found in the lookup table.')

    print(f"Extract weather files from archive")

    weather_archive = WEATHER_FILE_ARCHIVE(ldas)

    cmd = [
        SEVEN_ZIP,
        "e",
        f"{WEATHER_DIR}/{weather_archive}",
        f"-o{directory}",
    ]

    for w in df["Weather"].unique():
        cmd.append(w)

    cmd.append("-y")

    subprocess.run(
        cmd,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

    df.to_csv(f"{directory}/lut.csv", index=False)


    #    if len(region) != 3:
    #        cmd = f"{SEVEN_ZIP} a -sdel {region_strs}_{crop}/{ldas}_{crop}_{region_strs}_weather.7z ./{region_strs}_{crop}/*.weather"
    #        subprocess.run(
    #            cmd,
    #            shell=True,
    #            stdout=subprocess.DEVNULL,
    #            stderr=subprocess.DEVNULL,
    #        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter lookup tables using countries/regions")
    parser.add_argument(
        "--crop",
        required=True,
        choices=list(CROPS.keys()),
        help="Crop name",
    )
    parser.add_argument(
        "--ldas",
        required=True,
        choices=["GLDAS", "NLDAS", "gridMET"],
        help="LDAS type for weather files.",
    )
    parser.add_argument(
        "--directory",
        default="./",
        help="Directory to save output.",
    )
    args = parser.parse_args()

    main(args)
