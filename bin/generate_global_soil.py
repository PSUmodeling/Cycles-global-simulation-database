#!usr/bin/env python3
import argparse
import csv
import numpy as np
import os
import pandas as pd
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


def sample_soilgrids(rds, arrs, lats, lons, var, depth, lat, lon):
    val = list(rds[var + depth].sample([(lon, lat)]))[0][0]

    if val < 0:   # If the value is missing, use the nearest neighbor
        flag = True
        idx = ((np.abs(lats - lat)).argmin(), (np.abs(lons - lon)).argmin())
        radius = 1
        while radius <= 3:
            # Get the neighbors
            neighbors = arrs[var+depth][idx[0]-radius:idx[0]+radius+1, idx[1]-radius:idx[1]+radius+1]

            # If there is any valid neighbor, use it
            if np.any(neighbors[neighbors > 0]):
                val = np.mean(neighbors[neighbors > 0])
                break
            else:
                radius += 1
    else:
        flag = False

    val = val / SG_VARS[var] if val >= 0 else -999

    return val, flag


def main(crop):

    os.makedirs(SOIL_DIR, exist_ok=True)

    rds = {}    # Datasets for direct sampling
    arrs = {}   # Numpy arrays for gap filling
    for var in SG_VARS:
        for d in SG_DEPTHS:
            rds[var+d] = rasterio.open(GEOTIFFS["SoilGrids"](var, d))
            arrs[var+d] = rds[var+d].read(1)

            # Get the lats and lons of SoilGrids
            if 'lats' not in locals():
                height, width = arrs[var+d].shape
                cols, rows = np.meshgrid(np.arange(width), np.arange(height))
                xs, ys = rasterio.transform.xy(rds[var + d].transform, rows, cols)
                lons= np.array(xs)[0]
                lats = np.array([y[0] for y in np.array(ys)])

    df = pd.read_csv(
        HSG_SLOPE_CSV,
        index_col="GID",
    )

    for lu in LU_TYPES:
        print(crop, lu)

        os.makedirs(TEMP_DIR, exist_ok=True)

        # Read crop lookup table
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
                filled = []

                f.write(f"# Cycles soil file for {row['NAME_2'] + ', ' if row['NAME_2'] else ''}{row['NAME_1'] + ', ' if row['NAME_1'] else ''}{row['NAME_0']}\n#\n")
                f.write(f"# Clay, sand, organic matter and bulk density are obtained from SoilGrids.\n")
                f.write(f"# The data are sampled at Latitude {row['CropLat']}, Longitude {row['CropLon']}.\n")
                f.write(f"# NO3, NH4, and fractions of horizontal and vertical bypass flows are default empirical values.\n#\n")
                if hsg == -999:
                    f.write(f"# Hydrologic soil group MISSING DATA.\n")
                else:
                    f.write(f"# Hydrologic soil group {CURVE_NUMBERS[str(hsg)][1]}.\n")
                    f.write(f"# The curve number for row crops with straight row treatment is used, assuming good drainage. Please refer to the curve\n")
                    f.write(f"# number table in the README file if cover type, treatment, or drainage is different.\n\n")

                cn = CURVE_NUMBERS[str(hsg)][0] if hsg != -999 else -999
                f.write("%-15s\t%d\n" % ("CURVE_NUMBER", cn))

                f.write("%-15s\t%.2f\n" % ("SLOPE", slope))

                f.write("%-15s\t%s\n" % ("TOTAL_LAYERS", "6"))
                f.write(('%-7s\t'*12 + '%s\n') % (
                    "LAYER", "THICK", "CLAY", "SAND", "ORGANIC", "BD", "FC", "PWP", "SON", "NO3", "NH4", "BYP_H", "BYP_V"
                ))

                f.write(('%-7s\t'*12 + '%s\n') % (
                    "#", "m", "%", "%", "%", "Mg/m3", "m3/m3", "m3/m3", "kg/ha", "kg/ha", "kg/ha", "-", "-"
                ))

                for l in range(len(SG_DEPTHS)):
                    f.write("%-7d\t" % (l + 1))
                    f.write("%-7.2f\t" % SG_THICK[l])

                    clay, flag = sample_soilgrids(rds, arrs, lats, lons, "clay", SG_DEPTHS[l], lat, lon)
                    f.write("%-7d\t" % -999 if clay < 0 else "%-7.1f\t" % clay)
                    if flag:
                        filled.append(f'clay at {SG_DEPTHS[l]}')

                    sand, flag = sample_soilgrids(rds, arrs, lats, lons, "sand", SG_DEPTHS[l], lat, lon)
                    f.write("%-7d\t" % -999 if sand < 0 else "%-7.1f\t" % sand)
                    if flag:
                        filled.append(f'sand at {SG_DEPTHS[l]}')

                    om, flag = sample_soilgrids(rds, arrs, lats, lons, "soc", SG_DEPTHS[l], lat, lon)
                    f.write("%-7d\t" % -999 if om < 0 else "%-7.2f\t" % om)
                    if flag:
                        filled.append(f'organic matter at {SG_DEPTHS[l]}')

                    bd, flag = sample_soilgrids(rds, arrs, lats, lons, "bdod", SG_DEPTHS[l], lat, lon)
                    f.write("%-7d\t" % -999 if bd < 0 else "%-7.2f\t" % bd)
                    if flag:
                        filled.append(f'bulk density at {SG_DEPTHS[l]}')

                    f.write(('%-7d\t'*3 + '%-7.1f\t'*2 + '%-7.1f\t%.1f\n') % (
                        -999, -999, -999, DEFAULT_NO3[l], DEFAULT_NH4, 0.0, 0.0
                    ))

                if len(filled) > 0:
                    f.write("\n# The following values were missing and were filled in with values of adjacent grids:\n")
                    f.write(',\n'.join(['#   ' + x for x in filled]))
                    f.write('\n')

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
