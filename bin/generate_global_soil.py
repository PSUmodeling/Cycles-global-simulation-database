#!usr/bin/env python3
import argparse
import pandas as pd
import csv
import os
import rasterio
import subprocess
from setting import GEOTIFFS
from setting import CROPS
from setting import LU_TYPES
from setting import LOOKUP_CSV
from setting import SOIL_DIR
from setting import SOIL_ARCHIVE
from setting import SG_VARS
from setting import SG_DEPTHS
from setting import SG_THICK
from setting import DEFAULT_NO3
from setting import DEFAULT_NH4
from setting import TEMP_DIR
from setting import CURVE_NUMBERS
from setting import HSG_SLOPE_CSV
from setting import SEVEN_ZIP


def sample_soilgrids(rds, var, depth, lat, lon):
    val = list(rds[var + depth].sample([(lon, lat)]))[0][0]
    val /= SG_VARS[var]

    return val


def main(crop):

    os.makedirs(SOIL_DIR, exist_ok=True)

    rds = {}
    for var in SG_VARS:
        for d in SG_DEPTHS:
            rds[var+d] = rasterio.open(GEOTIFFS["SoilGrids"](var, d))

    df = pd.read_csv(
        HSG_SLOPE_CSV,
        index_col="GID",
    )

    for lu in LU_TYPES:
        # Read crop lookup table
        print(crop, lu)

        os.makedirs(TEMP_DIR, exist_ok=True)

        with open(LOOKUP_CSV(crop, lu)) as f:
            reader = csv.reader(f, delimiter=',')

            headers = next(reader)
            data = [{h:x for (h,x) in zip(headers,row)} for row in reader]

        for row in data:
            print(row["GID"], row["Soil"], row["CropLat"], row["CropLon"])
            lat = float(row["CropLat"])
            lon = float(row["CropLon"])
            hsg = int(df.loc[row["GID"]][f"{lu}HSG"])
            slope = df.loc[row["GID"]][f"{lu}Slope"]

            with open(f"{TEMP_DIR}/{row['Soil']}", "w") as f:

                f.write(f"# Cycles soil file for {row['NAME_2'] + ', ' if row['NAME_2'] else ''}{row['NAME_1'] + ', ' if row['NAME_1'] else ''}{row['NAME_0']}\n\n")
                f.write(f"# Soil physical parameters are sampled at Latitude {row['CropLat']}, Longitude {row['CropLon']}.\n")
                if hsg == -999:
                    f.write(f"# Hydrologic soil group MISSING DATA.\n")
                else:
                    f.write(f"# Hydrologic soil group {CURVE_NUMBERS[str(hsg)][1]}.\n")
                    f.write(f"# The curve number for row crops with straight row treatment is used, assuming good drainage. Please refer to the curve\n")
                    f.write(f"#   number table in the README file if cover type, treatment, or drainage is different.\n\n")

                cn = CURVE_NUMBERS[str(hsg)][0] if hsg != -999 else -999
                f.write("%-15s\t%d\n" % ("CURVE_NUMBER", cn))

                f.write("%-15s\t%.2f\n" % ("SLOPE", slope))

                f.write("%-15s\t%s\n" % ("TOTAL_LAYERS", "6"))
                f.write("%-7s\t" % "LAYER")
                f.write("%-7s\t" % "THICK")
                f.write("%-7s\t" % "CLAY")
                f.write("%-7s\t" % "SAND")
                f.write("%-7s\t" % "ORGANIC")
                f.write("%-7s\t" % "BD")
                f.write("%-7s\t" % "FC")
                f.write("%-7s\t" % "PWP")
                f.write("%-7s\t" % "SON")
                f.write("%-7s\t" % "NO3")
                f.write("%-7s\t" % "NH4")
                f.write("%-7s\t" % "BYP_H")
                f.write("%-7s\n" % "BYP_V")

                f.write("%-7s\t" % "#")
                f.write("%-7s\t" % "m")
                f.write("%-7s\t" % "%")
                f.write("%-7s\t" % "%")
                f.write("%-7s\t" % "%")
                f.write("%-7s\t" % "Mg/m3")
                f.write("%-7s\t" % "m3/m3")
                f.write("%-7s\t" % "m3/m3")
                f.write("%-7s\t" % "kg/ha")
                f.write("%-7s\t" % "kg/ha")
                f.write("%-7s\t" % "kg/ha")
                f.write("%-7s\t" % "-")
                f.write("%-7s\n" % "-")

                for l in range(len(SG_DEPTHS)):
                    f.write("%-7d\t" % (l + 1))
                    f.write("%-7.2f\t" % SG_THICK[l])

                    clay = sample_soilgrids(rds, "clay", SG_DEPTHS[l], lat, lon)
                    f.write("%-7d\t" % -999 if clay < 0 else "%-7.1f\t" % clay)

                    sand = sample_soilgrids(rds, "sand", SG_DEPTHS[l], lat, lon)
                    f.write("%-7d\t" % -999 if sand < 0 else "%-7.1f\t" % sand)

                    om = sample_soilgrids(rds, "soc", SG_DEPTHS[l], lat, lon)
                    f.write("%-7d\t" % -999 if om < 0 else "%-7.2f\t" % om)

                    bd = sample_soilgrids(rds, "bdod", SG_DEPTHS[l], lat, lon)
                    f.write("%-7d\t" % -999 if bd < 0 else "%-7.2f\t" % bd)

                    f.write("%-7d\t" % -999)
                    f.write("%-7d\t" % -999)
                    f.write("%-7d\t" % -999)

                    f.write("%-7.1f\t" % DEFAULT_NO3[l])
                    f.write("%-7.1f\t" % DEFAULT_NH4)

                    f.write("%-7s\t" % "0.0")
                    f.write("%-7s\n" % "0.0")

                f.write("\n")

        # Add to soil archive
        cmd = f"{SEVEN_ZIP} a -sdel {SOIL_ARCHIVE(crop, lu)} {TEMP_DIR}/*.soil"
        subprocess.run(
            cmd,
            shell=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )

        # Delete temp files
        cmd = [
            "rm",
            "-rf",
            TEMP_DIR,
        ]
        subprocess.run(
            cmd,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create soil files for different crops")
    parser.add_argument(
        "--crop",
        default="maize",
        choices=list(CROPS.keys()),
        help="Crop",
    )
    args = parser.parse_args()

    main(args.crop)
