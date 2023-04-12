#!/usr/bin/env python3

VERSION = "3.1"

CROPS = {
    "bean": "Pulses",
    "cassava": "Cassava",
    "lentil": "Pulses",
    "maize": "Maize",
    "millet": "Millet",
    "potato": "PotatoAndSweetpotato",
    "rice": "Rice",
    "sorghum": "Sorghum",
    "soybean": "Soybean",
    "sweetpotato": "PotatoAndSweetpotato",
    "wheat": "Wheat",
}

LU_TYPES = [
    "Rainfed",
    "Irrigated",
]

SG_VARS = {
    "bdod": 100.0,
    "clay": 10.0,
    "sand": 10.0,
    "soc": 58.0,
    #"cec": 10.0,
    #"cfvo": 10.0,
    #"nitrogen": 100000.0,
    #"ocd": 10.0,
    #"phh2o": 10.0,
    #"silt": 10.0,
}

SG_DEPTHS = [
    "0-5cm",
    "5-15cm",
    "15-30cm",
    "30-60cm",
    "60-100cm",
    "100-200cm",
]

SG_THICK = [
    0.05,
    0.1,
    0.15,
    0.3,
    0.4,
    1.0,
]

DEFAULT_NO3 = [
    10,
    10,
    7,
    4,
    4,
    4,
]

DEFAULT_NH4 = 1

CROP_VARS = [
    "HarvestedAreaFraction",
    "HarvestedAreaRainfed",
    "HarvestedAreaIrrigated",
    "Production",
    "ProductionRainfed",
    "ProductionIrrigated",
]

DATA_DIR = "./data"
LOOKUP_DIR = f"./crop_lookup_{VERSION}"
SOIL_DIR = f"./soil_{VERSION}"
TEMP_DIR = "./temp"

GADM_GPKG = f"{DATA_DIR}/gadm/gadm_410.gpkg"
GADM_SHP = f"{DATA_DIR}/gadm/gadm_410_l3.shp"
GADM_CSV = f"gadm_level3.csv"

HSG_SLOPE_CSV = f"hsg_slope.csv"

EARTHSTAT_CROP = lambda crop, var: f"{DATA_DIR}/HarvestedAreaYield175Crops_Geotiff/{crop}_HarvAreaYield_Geotiff/{crop}_{var}.tif"
GAEZ_CROP = lambda crop, type, var: f"{DATA_DIR}/GAEZ._2015_crop/GAEZAct2015_{var}_{crop}_{type}.tif"

GEOTIFFS = {
    "HarvestedAreaFraction": lambda crop: EARTHSTAT_CROP(crop, "HarvestedAreaFraction"),
    "HarvestedAreaRainfed": lambda crop: GAEZ_CROP(CROPS[crop], "Rainfed", "HarvArea"),
    "HarvestedAreaIrrigated": lambda crop: GAEZ_CROP(CROPS[crop], "Irrigated", "HarvArea"),
    "Production": lambda crop: EARTHSTAT_CROP(crop, "Production"),
    "ProductionRainfed": lambda crop: GAEZ_CROP(CROPS[crop], "Rainfed", "Production"),
    "ProductionIrrigated": lambda crop: GAEZ_CROP(CROPS[crop], "Irrigated", "Production"),
    "SoilGrids_org": lambda var, depth: f"{DATA_DIR}/SoilGrids/{var}/{var}_{depth}_mean_5000.tif",
    "SoilGrids": lambda var, depth: f"{DATA_DIR}/SoilGrids/{var}/{var}_{depth}_mean_5000_latlon.tif",
    "HSG": f"{DATA_DIR}/GlobalHSG/HYSOGs250m.tif",
    "slope": f"{DATA_DIR}/slope/slope_1KMmd_GMTEDmd.tif",
    "lc": f"{DATA_DIR}/ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7.tif",
}

FILTERED = f"gadm_filtered_{VERSION}.csv"

LOOKUP_CSV = lambda crop, lu: f"{LOOKUP_DIR}/{crop}_{lu.lower()}_global_lookup_{VERSION}.csv"
LOOKUP_CONUS_CSV = lambda crop, lu: f"{LOOKUP_DIR}/{crop}_{lu.lower()}_conus_lookup_{VERSION}.csv"

SOIL_ARCHIVE = lambda crop, lu: f"{SOIL_DIR}/{crop}_{lu.lower()}_global_soil_{VERSION}.7z"

GLDAS_MASK = f"{DATA_DIR}/GLDASp5_landmask_025d.nc4"
NLDAS_MASK = f"{DATA_DIR}/NLDAS_masks-veg-soil.nc4"

RAINFED_CROPLAND = [10, 11, 12]
IRRIGATED_CROPLAND = [20]
MOSAIC_CROPLAND = [30, 40]

CURVE_NUMBERS = {
    "1": [67, "A: low runoff potential (> 90% sand and < 10% clay)"],
    "2": [78, "B: moderately low runoff potential (50-90% sand and 10-20% clay)"],
    "3": [85, "C: moderately high runoff potential (< 50% sand and 20-40% clay)"],
    "4": [89, "D: high runoff potential (< 50% sand and > 40% clay)"],
    "11": [67,"A/D: high runoff potential unless drained (> 90% sand and < 10% clay)"],
    "12": [78,"B/D: high runoff potential unless drained (50-90% sand and 10-20% clay)"],
    "13": [85,"C/D: high runoff potential unless drained (< 50% sand and 20-40% clay)"],
    "14": [89,"D/D: high runoff potential unless drained (< 50% sand and > 40% clay)"],
}

SEVEN_ZIP = "7zzs"
COMPRESS = "./bin/compress.sh"
