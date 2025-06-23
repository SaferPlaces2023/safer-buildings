import os
import sys
import json
import click
import time
import logging
import requests
import datetime

import numpy as np
import pandas as pd

from concurrent.futures import ProcessPoolExecutor

from osgeo import gdal, ogr, osr
from shapely.wkt import loads
from shapely.geometry import box, Point, Polygon, MultiPolygon, LineString, MultiLineString
import geopandas as gpd

import leafmap

from . import _utils
from .module_log import Logger
from . import module_logo, module_args, module_retriever, module_flood, module_stats, module_s3, module_version


from dotenv import load_dotenv
load_dotenv()  


# TODO: --list-providers → shows detailed description of providers
# TODO: --summary → add metadata with aggregate stats ( n buildings, n flooded buildings, min,avg,max wd for each category (provider dependend) of buildings)


# DOC: Main function to compute flooded buildings

def compute_flood(
    water: str | None = None,
    building: str | None = None,
    wd_thresh: float = 0.5,
    bbox: tuple[float, float, float, float] | None = None,
    out: str | None = None,
    t_srs: str | None = None,
    provider: list[str] | None = None,
    filters: dict[str, list[dict[str, list|str|int]]] | None = None,
    only_flood: bool = False,
    stats: bool = False,
    summary: bool = False,
    
    # Additional parameters for CLI
    version: bool = False,
    debug: bool = False,
    verbose: bool = False
) -> str:
    
    """
    Main function to run the flooded buildings analysis.

    Parameters:
        waterdepth_filename (str): Path to the water depth raster file.
        buildings_filename (str | None): Path to the buildings vector file.
        bbox (tuple[float, float, float, float] | None): Bounding box for filtering.
        out (str | None): Output path for the results.
        t_srs (str | None): Target spatial reference system.
        providers (list[str] | None): List of data providers.
        feature_filters (dict[str, dict] | None): Filters for providers-features.

    Returns:
        str: Path to the output file containing flooded buildings.
    """
    
    out_cli = _utils.process_cli_args(
            version = version,
            debug = debug,
            verbose = verbose
        )
    if out_cli is not None:
        return out_cli
    

    # DOC: 0 — Print GECO logo
    module_logo.logo()
    
    
    # DOC: 1 — Validate args.
    Logger.debug("# Validating input arguments ...")
    validated_args = module_args.validate_args(
        waterdepth_filename = water,
        buildings_filename = building,
        wd_thresh = wd_thresh,
        bbox = bbox,
        out = out,
        t_srs = t_srs,
        provider = provider,
        feature_filters = filters,
        only_flood = only_flood,
        compute_stats = stats,
        compute_summary = summary
    )
    waterdepth_filename, buildings_filename, wd_thresh, bbox, out, t_srs, provider, feature_filters, only_flood, compute_stats, compute_summary = validated_args
    
    
    # DOC: 2 — Gather buildings
    Logger.debug(f'# Gather buildings data ...')
    provider_buildings = module_retriever.retrieve_buildings(
        buildings_filename=buildings_filename,
        bbox=bbox,
        provider=provider
    )
    
    
    # DOC: 3 — Polygonize waterdepth (compute one time and reuse for multiple providers)
    Logger.debug('# Processing water depth data (raster to polygon) ...')
    waterdepth_polygonized = _utils.polygonize_raster_valid_data(
        raster_filename=waterdepth_filename,
        band=1,
        mask_builder=lambda wd: wd > wd_thresh,    # Significant water depth threshold
        bbox=bbox
    )
    
    
    # DOC: 4 — Intersect buildings with water depth
    Logger.debug('# Intersecting buildings with water depth ...')
    flooded_buildings = module_flood.get_flooded_buildings(
        waterdepth_mask = waterdepth_polygonized,
        buildings = provider_buildings
    )
        
    
    # DOC: 5 — Filter features
    Logger.debug('# Filtering features ...')
    filtered_flooded_buildings = module_stats.filter_by_feature(
        gdf = flooded_buildings,
        feature_filters = feature_filters,
        only_flood = only_flood
    )
    Logger.debug(f"## Filtered {len(filtered_flooded_buildings)} buildings out from {len(flooded_buildings)}.")
    
    
    # DOC: 6 — Compute water depth stats over flooded buildings
    if compute_stats:
        Logger.debug('# Computing water depth stats over flooded buildings ...')
        filtered_flooded_buildings = module_stats.compute_wd_stats(
            waterdepth_filename=waterdepth_filename,
            waterdepth_mask=waterdepth_polygonized,
            buildings=filtered_flooded_buildings
        )
        Logger.debug("## Water depth stats computed for flooded buildings.")
        
    
    # DOC: 7 — Compute summary if requested
    if compute_summary:
        Logger.debug('# Computing summary statistics for flooded buildings ...')
        summary_stats = module_stats.compute_wd_summary(
            buildings=filtered_flooded_buildings,
            provider=provider,
            include_stats=compute_stats,
        )
        Logger.debug("## Summary statistics computed for flooded buildings.") 
    
        
    
    # DOC: 8 — Return results
    Logger.debug('# Preparing geojson output results ...')
    filtered_flooded_buildings = filtered_flooded_buildings.to_crs(t_srs)
    feature_collection = filtered_flooded_buildings.to_geo_dict()
    feature_collection['metadata'] = {
        'provider': provider,
        'buildings_count': len(filtered_flooded_buildings),
        'flooded_buildings_count': int(filtered_flooded_buildings['is_flooded'].sum()),
        'summary': summary_stats if compute_summary else None
    }
    feature_collection['crs'] = {
        "type": "name",
        "properties": {
            "name": f"urn:ogc:def:crs:{t_srs.replace(':', '::')}"  # REF: https://gist.github.com/sgillies/1233327 lines 256:271
        }
    }
    Logger.debug(f"## Buildings feature collection prepared with {len(feature_collection['features'])} features.")
        
    
    # DOC: 9 — Save results to file
    Logger.debug(f'# Saving results to {out} ...')
    if out.startswith('s3://'):
        out_tmp = _utils.temp_filename(ext='geojson', prefix='safer-buildings_out')
        with open(out_tmp, 'w') as f:
            json.dump(feature_collection, f, indent=2)
        module_s3.s3_upload(filename = out_tmp, uri = out)    
    else:
        with open(out, 'w') as f:
            json.dump(feature_collection, f, indent=2)
        
    
    Logger.debug(f"## Results saved to {out}")
    
    return feature_collection



# DOC: Main function to run the flooded buildings analysis from command line

@click.command()
@click.option('--water', type=str, required=False, help='Path to the water depth raster file.')
@click.option('--building','--buildings', type=str, default=None, help='Path to the buildings vector file.')
@click.option('--wd_thresh', type=float, default=0.5, help='Water depth threshold for significant flooding (default: 0.5).')
@click.option('--bbox', type=float, nargs=4, default=None, help='Bounding box (minx, miny, maxx, maxy). If None, the total bounds of the water depth raster will be used.')
@click.option('--out', type=str, default=None, help='Output path for the results.')
@click.option('--t_srs', type=str, default=None, help='Target spatial reference system (EPSG code). If None, CRS of water depth raster will be used.')
@click.option('--provider', type=str, default=None, help='Building data provider (one of OVERTURE, RER-REST/*).')
@click.option('--filters', type=str, default=None, help='Filters for providers-features in JSON format.')
@click.option('--only_flood', is_flag=True, required=False, default=False, help="Only return flooded buildings (default: False).")
@click.option('--stats', is_flag=True, required=False, default=False, help="Compute water depth statistics for flooded buildings.")
@click.option('--summary', is_flag=True, required=False, default=False, help="Returns an additional metadata field with aggregated statistic based on building type and class. If true, stats will be computed as well.")

@click.option("--version", is_flag=True, required=False, default=False,
              help="Print version.")
@click.option("--debug", is_flag=True, required=False, default=False,
              help="Debug mode.")
@click.option("--verbose", required=False, is_flag=True, type=click.BOOL, default=False,
              help="verbose mode.")
def main(
    water,
    building,
    wd_thresh,
    bbox,
    out,
    t_srs,
    provider,
    filters,
    only_flood,
    stats,
    summary,
    version,
    debug,
    verbose
):
    """
    Main function to run the flooded buildings analysis from command line.
    
    Examples:
    1. safer-buildings --water tests\rimini-wd.tif --provider OVERTURE --filters "[{'subtype':'education', 'class': ['kindergarten','school']}, {'class':'parking'}]" --only_flood --stats
    2. safer-buildings --water tests\rimini-wd.tif --provider RER-REST/28/31/40 --filters "[{'ORDINE_NORMALIZZATO': ['Scuola primaria', 'Nido d\'infanzia']}, {'ISTITUZIONE_SCOLASTICA_RIF': 'IC ALIGHIERI'}]" --summary
    3. safer-buildings --water s3://saferplaces.co/Directed/process_out/SaferBuildingsService/rimini-wd.tif --buildings s3://saferplaces.co/Directed/process_out/SaferBuildingsService/Data/buildings-default-area__rer-rest_overture.geojson --out s3://saferplaces.co/Directed/process_out/SaferBuildingsService/rimini-wd-buildings.geojson --provider RER-REST/1 --summary
    
    In first example the water depth file is 'tests/rimini-wd.tif', the OVERTURE provider is used, and buildings are filtered is (subtype in ['education'] AND class in ['kindergarten', 'school']) OR (class in ['parking']).
    """
    
    out_cli = _utils.process_cli_args(
        version = version,
        debug = debug,
        verbose = verbose
    )
    if out_cli is not None:
        sys.exit(0)
    
    Logger.debug("# Starting flooded buildings analysis ...")
    time_start = time.time()
    
    if filters:
        filters = json.loads(filters.replace("'", '"'))  # Convert single quotes to double quotes for JSON parsing
        
    Logger.debug("# CLI parameters:")
    Logger.debug(f"## Water depth file: {water}")
    Logger.debug(f"## Buildings file: {building}")
    Logger.debug(f"## Water depth threshold: {wd_thresh}")
    Logger.debug(f"## Bounding box: {bbox}")
    Logger.debug(f"## Output file: {out}")
    Logger.debug(f"## Target SRS: {t_srs}")
    Logger.debug(f"## Provider: {provider}")
    Logger.debug(f"## Feature filters: {filters}")
    Logger.debug(f"## Only flood: {only_flood}")
    Logger.debug(f"## Stats: {stats}")
    Logger.debug(f"## Summary: {summary}")
    
    result = compute_flood(
        water = water,
        building = building,
        wd_thresh = wd_thresh,
        bbox = tuple(bbox) if bbox else None,
        out = out,
        t_srs = t_srs,
        provider = provider,
        filters = filters,
        only_flood = only_flood,
        stats = stats,
        summary = summary,
        
        # Additional parameters for CLI
        version = version,
        debug = debug,
        verbose = verbose
    )
    
    time_end = time.time()
    
    if result is not None:
        Logger.debug(f"# Flooded buildings analysis completed in {time_end - time_start:.2f} seconds. Returned {len(result['features'])} features.")