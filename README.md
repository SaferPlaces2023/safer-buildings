# safer-buildings
Little procedures to get buildings at risk of flooding from a flood map.

## Installation option

##### For developers:
- `pip install .` no gdal install
- `pip install .[with_gdal_py37]` install gdal for python 3.7
- `pip install .[with_gdal_py39]` install gdal for python 3.9
- `pip install .[with_gdal_py310]` install gdal for python 3.10
- `pip install .[with_gdal_py311]` install gdal for python 3.11
- `pip install .[with_gdal_py312]` install gdal for python 3.12

##### By CLI:
- `pip install "git+https://github.com/SaferPlaces2023/safer-buildings.git#egg=safer-buildings[with_gdal_py***]"`

##### As dependency in your pyproj.toml-dependencies or setup.py-install_requires:
- `safer-buildings[with_gdal_py***] @ git+https://github.com/SaferPlaces2023/safer-buildings.git`


## Examples:
1. `safer-buildings --water <waterdepth.tif> --provider OVERTURE --filters "[{'subtype':'education', 'class': ['kindergarten','school']}, {'class':'parking'}] --only_flood --stats"`
2. `safer-buildings --water <waterdepth.tif> --provider REGIONE-EMILIA-ROMAGNA-30 --filters "[{'ORDINE_NORMALIZZATO': ['Scuola primaria', 'Nido d\'infanzia']}, {'ISTITUZIONE_SCOLASTICA_RIF': 'IC ALIGHIERI'}]"`
3. `safer-buildings --water s3://saferplaces.co/Directed/process_out/SaferBuildingsService/rimini-wd.tif --buildings s3://saferplaces.co/Directed/process_out/SaferBuildingsService/Data/buildings-default-area__rer-rest_overture.geojson --out s3://saferplaces.co/Directed/process_out/SaferBuildingsService/rimini-wd-buildings.geojson --provider RER-REST/1 --summary`

In first example the water depth file is 'tests/rimini-wd.tif', the OVERTURE provider is used, and buildings are filtered is `(subtype in ['education'] AND class in ['kindergarten', 'school']) OR (class in ['parking'])`.

### Parameters doc: 

| **Argument**                      | **Description**                                                                                                                                         |
|-------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------|
| **`--water`**, `--waterdepth`, `--wd`, `--water_filename` `TEXT` | Path to the water depth raster file.                                                                                                   |
| **`--building`**, `--buildings`, `--buildings_filename` `TEXT`   | Path to the buildings vector file.                                                                                                     |
| **`--wd_thresh`**, `--thresh` `FLOAT`      | Water depth threshold for significant flooding (default: `0.5`).                                                                     |
| **`--bbox`** `FLOAT...`                    | Bounding box (`minx`, `miny`, `maxx`, `maxy`). If None, the total bounds of the water depth raster will be used.                     |
| **`--out`** `TEXT`                         | Output path for the results.                                                                                                           |
| **`--t_srs`** `TEXT`                       | Target spatial reference system (EPSG code). If None, CRS of building will be used if provided, otherwise CRS of water depth raster. |
| **`--provider`** `TEXT`                    | Building data provider (one of `OVERTURE`, `RER-REST/*`).                                                                              |
| **`--filters`** `TEXT`                     | Filters for provider features in JSON format.                                                                                         |
| **`--only_flood`**                       | Only return flooded buildings (default: `False`).                                                                                      |
| **`--stats`**                            | Compute water depth statistics for flooded buildings.                                                                                 |
| **`--summary`**                          | Add aggregated statistics metadata based on building type and class.                                               |
| **`--summary_on`** `TEXT,...`                         | Column(s) _(separated by commas – no spaces allowed)_ to compute summary statistics on. If `None`, summary will be computed on all flooded buildings. If not provided and provider is `OVERTURE`, `'subtype'` will be used, if provider is `RER-REST` then `'rer_class'` will be used.                                               |
| **`--out_geojson`**                      | Output results in GeoJSON format (default: `False`). If `True`, output will be a GeoJSON feature collection, otherwise a JSON with references. |
| **`--version`**, `-v`                    | Print version.                                                                                                                         |
| **`--debug`**                            | Enable debug mode.                                                                                                                     |
| **`--verbose`**                          | Enable verbose mode.                                                                                                                   |
| `--help`                             | Show this message and exit.                                                                                                            |

#### _Please note that_:
- **bold arguments** are used as input state in AWS Lambda. <br>
- `TYPE ...` arguments need to be expressed as `list[type]` in AWS Lambda.
