# safer-buildings
Little procedures to get buildings at risk of flooding from a flood map.

Examples:
1. `safer-buildings --wd <waterdepth.tif> --provider OVERTURE --filters "[{'subtype':'education', 'class': ['kindergarten','school']}, {'class':'parking'}] --only_flood --stats"`
2. `safer-buildings --wd <waterdepth.tif> --provider REGIONE-EMILIA-ROMAGNA-30 --filters "[{'ORDINE_NORMALIZZATO': ['Scuola primaria', 'Nido d\'infanzia']}, {'ISTITUZIONE_SCOLASTICA_RIF': 'IC ALIGHIERI'}]"`

In first example the water depth file is 'tests/rimini-wd.tif', the OVERTURE provider is used, and buildings are filtered is `(subtype in ['education'] AND class in ['kindergarten', 'school']) OR (class in ['parking'])`.
