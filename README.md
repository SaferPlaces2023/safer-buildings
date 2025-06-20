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
