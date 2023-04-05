# Cycles global simulation database version 3.1

The Cycles global simulation database version 3.1 provides crop lookup tables, Cycles soil files, and Cycles weather files (1979 to present in CONUS, 2000 to present globally).
This version supports the simulations of major crops in any Level-3 (e.g., county level) administrative region in the world.
The crops that can be simulated include **bean, cassava, lentil, maize, millet, potato, rice, sorghum, soybean, sweet potato, and wheat**.
Each crop is classified into two categories based on irrigation types: rainfed and irrigated.
All versions prior to v2.0 were developed by Dr. Lorne Leonard.

## Data Sources

- [GAEZ+_2015 Crop Harvest Area](https://doi.org/10.7910/DVN/KAGRFI)
- [GAEZ+_2015 Crop Production](https://doi.org/10.7910/DVN/KJFUO1)
- [Harvested Area and Yield for 175 Crops](http://www.earthstat.org/harvested-area-yield-175-crops/)
- [GADM v4.1](https://gadm.org/data.html)
- [ESA/CCI Land Cover Maps - v2.0.7](http://maps.elie.ucl.ac.be/CCI/viewer/download.php#ftp_dwl)
- [SoilGrids aggregated 5000-m](https://files.isric.org/soilgrids/latest/data_aggregated/5000m/)
- [Global Hydrologic Soil Groups (HYSOGs250m) for Curve Number-Based Runoff Modeling](https://doi.org/10.3334/ORNLDAAC/1566)
- [Curve number tables](https://www.hec.usace.army.mil/confluence/hmsdocs/hmstrm/cn-tables)
- [GMTED2010 1-km median slope](https://www.earthenv.org/topography)

## Crop lookup table files

The `crop_lookup_3.1` directory contains crop lookup tables that include all 3rd-level (or above) administrative regions where major crops are harvested, along with the names of corresponding weather files and soil files for the regions.
The lookup tables are provided in `csv` format, and are named using the convention
`[crop name]_[irrigation type]_[range]_lookup_[file version].csv`.
For example, `maize_irrigated_global_lookup_3.1.csv` can be interpreted as follows:
- [crop name] = maize
- [irrigation type] = irrigated
- [range] = global
- [File version] = 3.1

Each lookup table file is structured as:
Column          | Description
--------------- | -----------------------------------------------------------
GID             | Unique id representing Administrative region of the world
NAME_0          | Country
NAME_1          | Level 1 sub-division (e.g., state)
NAME_2          | Level 2 sub-division (e.g., county)
Lat             | Latitude of region centroid (degree)
Lon             | Longitude of region centroid (degree)
AreaKm2         | Region area in km<sup>2</sup>
AreaFraction    | Harvested area fraction of the crop
CropLat         | Latitude of maximum crop production in the area
CropLon         | Longitude of maximum crop production in the area
Weather         | Weather file name
Soil            | Soil file name

Global look-up tables map regions to GLDAS grids, and CONUS look-up tables match them to NLDAS grids.

## Soil file archives

Soil file archives provide global soil files that support the simulations of major crops, and can be found in the `soil_3.1` directory.
These files are zipped and follow the naming convention `[crop name]_[irrigation type]_global_soil_[file version].7z`.

The soil files contain comments indicating the hydrologic soil group to which the soil belongs.
The curve number for row crops with straight row treatment is used, assuming good drainage.
However, if the cover type, treatment, or drainage is different, please refer to the table below to find the appropriate curve numbers.

Cover type                                           | Treatment  | Hydrologic condition | Group A | Group B | Group C | Group D
-----------------------------------------------------|------------|----------------------|---------|---------|---------|--------
Fallow                                               | Bare soil  |                      |   77    |   86    |   91    |   94
Fallow                                               | CR         | Poor                 |   76    |   85    |   90    |   93
Fallow                                               | CR         | Good                 |   74    |   83    |   88    |   90
Row crops                                            | SR         | Poor                 |   72    |   81    |   88    |   91
Row crops                                            | SR         | Good                 |   67    |   78    |   85    |   89
Row crops                                            | SR + CR    | Poor                 |   71    |   80    |   87    |   90
Row crops                                            | SR + CR    | Good                 |   64    |   75    |   82    |   85
Row crops                                            | C          | Poor                 |   70    |   79    |   84    |   88
Row crops                                            | C          | Good                 |   65    |   75    |   82    |   86
Row crops                                            | C + CR     | Poor                 |   69    |   78    |   83    |   87
Row crops                                            | C + CR     | Good                 |   64    |   74    |   81    |   85
Row crops                                            | C & T      | Poor                 |   66    |   74    |   80    |   82
Row crops                                            | C & T      | Good                 |   62    |   71    |   78    |   81
Row crops                                            | C & T + CR | Poor                 |   65    |   73    |   79    |   81
Row crops                                            | C & T + CR | Good                 |   61    |   70    |   77    |   80
Small grain                                          | SR         | Poor                 |   65    |   76    |   84    |   88
Small grain                                          | SR         | Good                 |   63    |   75    |   83    |   87
Small grain                                          | SR + CR    | Poor                 |   64    |   75    |   83    |   86
Small grain                                          | SR + CR    | Good                 |   60    |   72    |   80    |   84
Small grain                                          | C          | Poor                 |   63    |   74    |   82    |   85
Small grain                                          | C          | Good                 |   61    |   73    |   81    |   84
Small grain                                          | C + CR     | Poor                 |   62    |   73    |   81    |   84
Small grain                                          | C + CR     | Good                 |   60    |   72    |   80    |   83
Small grain                                          | C & T      | Poor                 |   61    |   72    |   79    |   82
Small grain                                          | C & T      | Good                 |   59    |   70    |   78    |   81
Small grain                                          | C & T + CR | Poor                 |   60    |   71    |   78    |   81
Small grain                                          | C & T + CR | Good                 |   58    |   69    |   77    |   80
Close-seeded or broadcast legumes or rotation meadow | SR         | Poor                 |   66    |   77    |   85    |   89
Close-seeded or broadcast legumes or rotation meadow | SR         | Good                 |   58    |   72    |   81    |   85
Close-seeded or broadcast legumes or rotation meadow | C          | Poor                 |   64    |   75    |   83    |   85
Close-seeded or broadcast legumes or rotation meadow | C          | Good                 |   55    |   69    |   78    |   83
Close-seeded or broadcast legumes or rotation meadow | C & T      | Poor                 |   63    |   73    |   80    |   83
Close-seeded or broadcast legumes or rotation meadow | C & T      | Good                 |   51    |   67    |   76    |   80

**Note:**

**CR: Crop residue cover**

**SR: Straight row**

**C: Contoured**

**T: Terraced**

## Weather file archives

Weather file archives are stored in the `weather` directory.
The `NLDAS_CONUS_1979-2022.7z` archive contains 52,476 Cycles weather files for the CONUS region, generated from the primary forcing data for Phase 2 of the North American Land Data Assimilation System (NLDAS-2).
The `GLDAS_2000-2022.7z` archive contains all 47221 Cycles weather files that appear in the global look-up tables, generated from the primary forcing data for the Global Land Data Assimilation System (GLDAS).
The weather files follow the naming convention `[LDAS]_[lat][N or S]_[lon][E or W].weather`, where `[LDAS]` can be either GLDAS or NLDAS, and `[lat]` and `[lon]` refer to the latitude and longitude of the corresponding grids.
