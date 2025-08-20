import os
import sys
import json
import click
import time
import pprint
import logging
import requests
import datetime
import traceback

import numpy as np
import pandas as pd

from concurrent.futures import ProcessPoolExecutor

from osgeo import gdal, ogr, osr
from shapely.wkt import loads
from shapely.geometry import box, Point, Polygon, MultiPolygon, LineString, MultiLineString
import geopandas as gpd

import leafmap
 
from . import _consts, _utils, module_s3, module_version
from . import module_logo, module_args, module_retriever, module_flood, module_stats, module_add_ops, module_outputs
from .module_args import _ARG_NAMES
from .module_log import Logger, is_debug_mode
from .module_exception import CustomException, ArgsException, RetrieverException, FloodException, StatsException, AddOpsException, OutputsException

from dotenv import load_dotenv
load_dotenv()  


# TODO: --list-providers → shows detailed description of providers


# DOC: Main function to compute flooded buildings

def compute_flood(
    water: str | None = None,
    building: str | None = None,
    wd_thresh: float = 0.5,
    bbox: tuple[float, float, float, float] | None = None,
    out: str | None = None,
    t_srs: str | None = None,
    provider: list[str] | None = None,
    flood_mode: str | None = None,
    filters: dict[str, list[dict[str, list|str|int]]] | None = None,
    only_flood: bool = False,
    stats: bool = False,
    summary: bool = False,
    summary_on: str | list[str] | None = None,
    add_ops: list[str] | None = None,
    out_geojson: bool = False,
    
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
    
    try:
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
        try:
            Logger.debug("# Validating input arguments ...")
            validated_args = module_args.validate_args(
                waterdepth_filename = water,
                buildings_filename = building,
                wd_thresh = wd_thresh,
                bbox = bbox,
                out = out,
                t_srs = t_srs,
                provider = provider,
                flood_mode = flood_mode,
                feature_filters = filters,
                only_flood = only_flood,
                compute_stats = stats,
                compute_summary = summary,
                summary_on = summary_on,
                add_ops = add_ops,
                out_geojson = out_geojson,
            )
            waterdepth_filename, buildings_filename, wd_thresh, bbox, out, t_srs, provider, flood_mode, feature_filters, only_flood, compute_stats, compute_summary, summary_on, add_ops, out_geojson = validated_args
        except Exception as e:
            raise ArgsException.from_exception(e)
        
        
        # DOC: 2 — Gather buildings
        try:
            Logger.debug(f'# Gather buildings data ...')
            provider_buildings = module_retriever.retrieve_buildings(
                buildings_filename=buildings_filename,
                bbox=bbox,
                provider=provider
            )
        except Exception as e:
            raise RetrieverException.from_exception(e)
        
        
        # DOC: 3 — Polygonize waterdepth (compute one time and reuse for multiple providers)
        try:
            Logger.debug('# Processing water depth data (raster to polygon) ...')
            waterdepth_polygonized = module_flood.poligonyze_waterdepth(
                waterdepth_filename = waterdepth_filename,
                wd_threshold = wd_thresh,
                bbox = bbox
            )
        except Exception as e:
            raise FloodException.from_exception(e)
        
        
        # DOC: 4 — Intersect buildings with water depth
        try:
            Logger.debug('# Intersecting buildings with water depth ...')
            flooded_buildings = module_flood.get_flooded_buildings(
                waterdepth_gdf = waterdepth_polygonized,
                buildings = provider_buildings,
                flood_mode = flood_mode
            )
        except Exception as e:
            raise FloodException.from_exception(e)
            
        
        # DOC: 5 — Filter features
        try:
            Logger.debug('# Filtering features ...')
            flooded_buildings = module_stats.filter_by_feature(
                gdf = flooded_buildings,
                feature_filters = feature_filters,
                only_flood = only_flood
            )
            Logger.debug(f"## Filtered {len(flooded_buildings)} buildings out from {len(flooded_buildings)}.")
        except Exception as e:
            raise StatsException.from_exception(e)
        
        # DOC: 6 — Compute water depth stats over flooded buildings
        try:
            if compute_stats:
                Logger.debug('# Computing water depth stats over flooded buildings ...')
                flooded_buildings = module_stats.compute_wd_stats(
                    waterdepth_filename=waterdepth_filename,
                    waterdepth_mask=waterdepth_polygonized,
                    waterdepth_thresh=wd_thresh,
                    buildings=flooded_buildings,
                    flood_mode=flood_mode
                )
                Logger.debug("## Water depth stats computed for flooded buildings.")
        except Exception as e:
            raise StatsException.from_exception(e)
            
        
        # DOC: 7 — Compute summary if requested
        try:
            summary_stats_ouput_data = dict()
            if compute_summary:
                Logger.debug('# Computing summary statistics for flooded buildings ...')
                summary_stats = module_stats.compute_wd_summary(
                    buildings=flooded_buildings,
                    summary_on=summary_on,
                    provider=provider,
                    include_stats=compute_stats,
                )
                summary_stats_ouput_data = { 'summary': summary_stats }
                Logger.debug("## Summary statistics computed for flooded buildings.") 
        except Exception as e:
            raise StatsException.from_exception(e)
        
        # DOC: 7.1 — Drop other geometry columns
        flooded_buildings = flooded_buildings.drop(columns = set(flooded_buildings.columns) & set([_consts._COL_FLOOD_ROI, _consts._COL_FLOOD_AREA, _consts._COL_FLOOD_GEOMETRY, _consts._COL_FLOOD_COORDS]))

        # DOC: 8 — Run additional operations if any
        try:
            add_ops_output_data = dict()
            if add_ops is not None:
                for op in add_ops:
                    Logger.debug(f'# Running additional operation: {op.name}...')
                    # DOC: Handle Additional Operation execution with its proper arguments
                    if op.name == module_add_ops.NearbyPumps.name:
                        nearby_pumps_output = op( ** {
                            'gdf_buildings': flooded_buildings,
                            'gdf_water_depth': waterdepth_polygonized,
                            'bbox': bbox,
                            't_srs': t_srs,
                        })
                        flooded_buildings, nearby_pumps_collection = nearby_pumps_output
                        add_ops_output_data[op.name] = nearby_pumps_collection
                    elif op.name == module_add_ops.AlertMethod.name:
                        alert_method_output = op( ** {
                            'gdf_buildings': flooded_buildings,
                            'gdf_water_depth': waterdepth_polygonized,
                            'bbox': bbox,
                            't_srs': t_srs,
                        })
                        flooded_buildings, alert_method_collection = alert_method_output
                        add_ops_output_data[op.name] = alert_method_collection
                    elif op.name == module_add_ops.GatesGuard.name:
                        gates_guard_output = op( ** {
                            'gdf_buildings': flooded_buildings,
                            'gdf_water_depth': waterdepth_polygonized,
                            'bbox': bbox,
                            't_srs': t_srs,
                        })
                        flooded_buildings, gates_collection = gates_guard_output
                        # add_ops_output_data[op.name] = gates_collection
        except Exception as e:
            raise AddOpsException.from_exception(e)


        # DOC: 9 — Return results
        try:
            Logger.debug('# Preparing geojson output results ...')
            feature_collection = module_outputs.prepare_feture_collection(
                flooded_buildings = flooded_buildings,
                t_srs = t_srs,
                provider = provider,
                summary_stats = summary_stats_ouput_data,
                add_ops_output_data = add_ops_output_data
            )   
        except Exception as e:
            raise OutputsException.from_exception(e)     
            
        
        # DOC: 10 — Save results to file
        try:
            Logger.debug(f'# Saving results to {out} ...')
            module_outputs.save_results(
                feature_collection = feature_collection,
                out = out
            )
            for add_op_name, add_op_collection in add_ops_output_data.items():
                Logger.debug(f'## Saving additional operation output: {add_op_name} ...')
                module_outputs.save_results(
                    feature_collection = add_op_collection,
                    out = out,
                    out_postfix = add_op_name
                )
        except Exception as e:
            raise OutputsException.from_exception(e)
        
        
        # DOC: 11 — Return output
        try: 
            Logger.debug('# Returning output ...')
            output = module_outputs.prepare_output(
                feature_collection = feature_collection,
                out = out,
                out_geojson = out_geojson,
                compute_summary = compute_summary
            )
        except Exception as e:
            raise OutputsException.from_exception(e)
        
        return output

    except CustomException as e:
        exception_info = e.to_dict()
        Logger.error(f"\nException occurred: \n{pprint.pformat(exception_info)}\n")
        return exception_info
    
    except Exception as e:
        custom_ex = CustomException.from_exception(e)
        custom_ex.message = f"An error occurred during the flooded buildings analysis: {custom_ex.message}"
        exception_info = custom_ex.to_dict()
        Logger.error(f"\nException occurred: \n{pprint.pformat(exception_info)}\n")
        return exception_info
    
    finally:
        # DOC: 12 — Clean up temporary files created in this run istance
        Logger.debug("# Cleaning up temporary files ...")
        _utils.clean_temp_files(from_garbage_collection=True)


# DOC: Main function to run the flooded buildings analysis from command line

@click.command()
@click.option(
    *_ARG_NAMES.WATER,
    type=str, required=False, help='Path to the water depth raster file.'
)
@click.option(
    *_ARG_NAMES.BUILDINGS,
    type=str, default=None, help='Path to the buildings vector file.'
)
@click.option(
    *_ARG_NAMES.WD_THRESH,
    type=float, default=0.5, help='Water depth threshold for significant flooding (default: 0.5).'
)
@click.option(
    *_ARG_NAMES.BBOX,
    callback=lambda ctx, param, value: tuple(value) if value else None,
    type=float, nargs=4, default=None, help='Bounding box (minx, miny, maxx, maxy). If None, the total bounds of the water depth raster will be used.'
)
@click.option(
    *_ARG_NAMES.OUT,
    type=str, default=None, help='Output path for the results.'
)
@click.option(
    *_ARG_NAMES.T_SRS,
    type=str, default=None, help='Target spatial reference system (EPSG code). If None, CRS of building will be used if file is provided otherwise CRS of water depth raster will be used.'
)
@click.option(
    *_ARG_NAMES.PROVIDER,
    type=str, default=None, help='Building data provider (one of OVERTURE, RER-REST/*, VENEZIA-WFS/*, VENEZIA-WFS-CRITICAL-SITES).'
)
@click.option(
    *_ARG_NAMES.FLOOD_MODE,
    type=str, default=None, help='Where to search for flooded buildings/areas. Valid values are \'buffer\', \'in-area\', \'all\'. When \'buffer\' use buffer look for flood around buildings, when \'in-area\' look for flood inside geometry, when \'all\' look for flood in both ways. If None, the default value is \'buffer\'.'
)
@click.option(
    *_ARG_NAMES.FILTERS,
    callback = lambda ctx, param, value: json.loads(value.replace("'", '"')) if value else None,
    type=str, default=None, help='Filters for providers-features in JSON format.'
)
@click.option(
    *_ARG_NAMES.ONLY_FLOOD,
    is_flag=True, required=False, default=False, help="Only return flooded buildings (default: False)."
)
@click.option(
    *_ARG_NAMES.STATS,
    is_flag=True, required=False, default=False, help="Compute water depth statistics for flooded buildings."
)
@click.option(
    *_ARG_NAMES.SUMMARY,
    is_flag=True, required=False, default=False, help="Returns an additional metadata field with aggregated statistic based on building type and class. If true, stats will be computed as well."
)
@click.option(
    *_ARG_NAMES.SUMMARY_ON,
    callback = lambda ctx, param, value: value.split(',') if value is not None else None,
    required=False, default=None, help="Column(s) (separated by commas — no spaces allowed) to compute summary statistics on. If None, summary will be computed on all flooded buildings. If not provided and provider is OVERTURE, 'subtype' will be used, if provider is RER-REST then 'service_class' will be used."
)
@click.option(
    *_ARG_NAMES.ADD_OPS,
    callback = lambda ctx, param, value: _ARG_NAMES.parse_add_ops(value),
    required=False, default=None, help=f"Additional operations to perform on the results. Implemented operations {' '.join(list(module_add_ops._ADD_OPS.keys()))}, but avaliable ones depend on selected provider and geographical area."
)
@click.option(
    *_ARG_NAMES.OUT_GEOJSON,
    is_flag=True, required=False, default=False, help="Output results in GeoJSON format (default: False). By default is setted to False due to large dimensions. If true, output will be a GeoJSON feature collection, if False, it will be a json with GeoJSON reference and metadata."
)

@click.option("--version", "-v", is_flag=True, required=False, default=False, help="Print version.")
@click.option("--debug", is_flag=True, required=False, default=False, help="Debug mode.")
@click.option("--verbose", required=False, is_flag=True, type=click.BOOL, default=False, help="verbose mode.")
def main(
    water,
    building,
    wd_thresh,
    bbox,
    out,
    t_srs,
    provider,
    flood_mode,
    filters,
    only_flood,
    stats,
    summary,
    summary_on,
    add_ops,
    out_geojson,
    
    version,
    debug,
    verbose
):
    """
    Main function to run the flooded buildings analysis from command line.
    
    Examples:
    1. safer-buildings --water tests\rimini-wd.tif --provider OVERTURE --filters "[{'subtype':'education', 'class': ['kindergarten','school']}, {'class':'parking'}]" --only_flood --stats
    2. safer-buildings --water tests\rimini-wd.tif --provider RER-REST/28/31/40 --filters "[{'ORDINE_NORMALIZZATO': ['Scuola primaria', 'Nido d\'infanzia']}, {'ISTITUZIONE_SCOLASTICA_RIF': 'IC ALIGHIERI'}]" --summary
    3. safer-buildings --water s3://saferplaces.co/Directed/process_out/SaferBuildingsService/rimini-wd.tif --buildings s3://saferplaces.co/Directed/process_out/SaferBuildingsService/Data/buildings-default-area__rer-rest_overture.geojson --out s3://saferplaces.co/Directed/process_out/SaferBuildingsService/rimini-wd-buildings.geojson --provider RER-REST --summary
    4. safer-buildings --water s3://saferplaces.co/api_data/vimercate/waterdepths/rain-100/water_rain-100.tif --buildings s3://saferplaces.co/api_data/vimercate/building.shp --out s3://saferplaces.co/api_data/vimercate/flooded-building.geojson --provider OVERTURE --summary
    
    In first example the water depth file is 'tests/rimini-wd.tif', the OVERTURE provider is used, and buildings are filtered is (subtype in ['education'] AND class in ['kindergarten', 'school']) OR (class in ['parking']).
    """
    
    out_cli = _utils.process_cli_args(
        version = version,
        debug = debug,
        verbose = verbose
    )
    if out_cli is not None:
        sys.exit(0)
    
    Logger.debug(f"# Starting flooded buildings analysis ... PID: {os.getpid()}")
    time_start = time.time()
    
    Logger.debug("# CLI parameters:")
    Logger.debug(f"## Water depth file: {water}")
    Logger.debug(f"## Buildings file: {building}")
    Logger.debug(f"## Water depth threshold: {wd_thresh}")
    Logger.debug(f"## Bounding box: {bbox}")
    Logger.debug(f"## Output file: {out}")
    Logger.debug(f"## Target SRS: {t_srs}")
    Logger.debug(f"## Provider: {provider}")
    Logger.debug(f"## Flood mode: {flood_mode}")
    Logger.debug(f"## Feature filters: {filters}")
    Logger.debug(f"## Only flood: {only_flood}")
    Logger.debug(f"## Stats: {stats}")
    Logger.debug(f"## Summary: {summary}")
    Logger.debug(f"## Summary on: {summary_on}")
    Logger.debug(f"## Additional operations: {add_ops}")
    Logger.debug(f'## Output GeoJSON: {out_geojson}')
    
    result = compute_flood(
        water = water,
        building = building,
        wd_thresh = wd_thresh,
        bbox = bbox,
        out = out,
        t_srs = t_srs,
        provider = provider,
        flood_mode = flood_mode,
        filters = filters,
        only_flood = only_flood,
        stats = stats,
        summary = summary,
        summary_on = summary_on,
        add_ops = add_ops,
        out_geojson = out_geojson,
        
        # Additional parameters for CLI
        version = version,
        debug = debug,
        verbose = verbose
    )
    
    time_end = time.time()
    
    if result is not None:
        Logger.debug(f"# Flooded buildings analysis completed in {time_end - time_start:.2f} seconds.")
        
    return result