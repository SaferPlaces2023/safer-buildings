import os
import re
import datetime
from typing import Any

from osgeo import osr

from shapely.geometry import box
import geopandas as gpd

from . import _utils, _consts, filesystem
from . import module_add_ops, module_s3

from .module_log import Logger



class _ARG_NAMES():
    WATER = ['--water', '--waterdepth', '--wd', '--water_filename']
    BUILDINGS = ['--building', '--buildings', '--buildings_filename']
    WD_THRESH = ['--wd_thresh', '--thresh']
    BBOX = ['--bbox']
    OUT = ['--out']
    T_SRS = ['--t_srs']
    PROVIDER = ['--provider']
    FLOOD_MODE = ['--flood_mode']
    FILTERS = ['--filters']
    ONLY_FLOOD = ['--only_flood']
    STATS = ['--stats']
    SUMMARY = ['--summary']
    SUMMARY_ON = ['--summary_on']
    ADD_OPS = ['--add_ops']
    OUT_GEOJSON = ['--out_geojson']


    def parse_add_ops(value: str) -> dict[str, dict[str, Any]] | None:

        def parse_op(op_str):
            # DOC: Pattern to match the operation name and attributes
            pattern = r'^(\w+)(?:\[(.*?)\])?$'
            attr_pattern = r'(\w+)\s*=\s*(?:"([^"]*)"|([^;]*))'

            match = re.match(pattern, op_str)
            if not match:
                return {}

            key = match.group(1)
            attr_str = match.group(2)

            attributes = {}
            if attr_str:
                for attr_match in re.finditer(attr_pattern, attr_str):
                    attr_name = attr_match.group(1)
                    attr_value = attr_match.group(2) or attr_match.group(3)
                    attributes[attr_name] = attr_value

            return {key: attributes}
        
        if value is None:
            return None
        ops_str = value.split('|')
        ops_dict = {}
        for op_str in ops_str:
            ops_dict.update(parse_op(op_str.strip()))
        return ops_dict



def validate_args(
    waterdepth_filename: str,
    buildings_filename: str | None = None,
    wd_thresh: float = 0.5,
    bbox: tuple[float, float, float, float] | None = None,
    out: str | None = None,
    t_srs: str | None = None,
    provider: list[str] | None = None,
    flood_mode: str | None = None,
    feature_filters: dict[str, dict] | None = None,
    only_flood: bool = False,
    compute_stats: bool = False,
    compute_summary: bool = False,
    summary_on: str | list[str] | None = None,
    add_ops: list[str] | None = None,
    out_geojson: bool = False
) -> tuple[str, str | None, float, gpd.GeoDataFrame, str, str, str, list[dict[str, list]], bool, bool, bool]:
    
    """
    Validate the input arguments for the flooded buildings analysis.
    """
        
    if waterdepth_filename is None:
        raise ValueError(f"water must be provided. Check water argument ({_ARG_NAMES.WATER})")
    if type(waterdepth_filename) is not str:
        raise TypeError(f"water must be a string. Check water argument ({_ARG_NAMES.WATER})")
    if not waterdepth_filename.startswith('s3://') and os.path.isfile(waterdepth_filename) is False:
        raise FileNotFoundError(f"Water depth file not found: {waterdepth_filename}. Check water argument ({_ARG_NAMES.WATER})")
    if waterdepth_filename.startswith('s3://'):
        waterdepth_filename_tmp = _utils.temp_filename(ext = _utils.justext(waterdepth_filename), prefix='safer-buildings_waterdepth')
        dl_status = module_s3.s3_download(uri = waterdepth_filename, fileout = waterdepth_filename_tmp)
        if dl_status is None:
            raise FileNotFoundError(f"Failed to download water depth file from S3: {waterdepth_filename}. Check water argument ({_ARG_NAMES.WATER})")
        waterdepth_filename = waterdepth_filename_tmp

    _consts.init_epsg_utmxx(_utils.epsg_utm_from_raster(waterdepth_filename))
    Logger.debug(f"## EPSG UTMxx initialized from waterdepth: {_consts._EPSG_UTMxx}")
    
    if buildings_filename is not None:
        if type(buildings_filename) is not str:
            raise TypeError(f"buildings must be a string. Check buildings arguments ({_ARG_NAMES.BUILDINGS})")
        if not buildings_filename.startswith('s3://') and os.path.isfile(buildings_filename) is False:
            raise FileNotFoundError(f"Buildings file not found: {buildings_filename}. Check buildings arguments ({_ARG_NAMES.BUILDINGS})")
        if buildings_filename.startswith('s3://'):
            buildings_filename_tmp = _utils.temp_filename(ext = _utils.justext(buildings_filename), prefix='safer-buildings_buildings')
            buildings_filename_uri_aux = filesystem.get_aux_files(buildings_filename)
            buildings_filename_tmp_aux = filesystem.get_aux_files(buildings_filename_tmp)
            for bf_uri_aux, bf_tmp_aux in zip(buildings_filename_uri_aux, buildings_filename_tmp_aux):
                dl_status = module_s3.s3_download(uri = bf_uri_aux, fileout=bf_tmp_aux)
                if dl_status is None:
                    raise FileNotFoundError(f"Failed to download buildings auxiliary files from S3: {bf_uri_aux}. Check buildings arguments ({_ARG_NAMES.BUILDINGS})")
            dl_status = module_s3.s3_download(uri = buildings_filename, fileout = buildings_filename_tmp)
            if dl_status is None:
                raise FileNotFoundError(f"Failed to download buildings file from S3: {buildings_filename}. Check buildings arguments ({_ARG_NAMES.BUILDINGS})")
            buildings_filename = buildings_filename_tmp
        
    if wd_thresh is None:
        wd_thresh = 0.5
    if type(wd_thresh) not in (int, float):
        raise TypeError(f"wd_thresh must be a float or int. Check water depth threshold argument ({_ARG_NAMES.WD_THRESH})")
    if wd_thresh < 0:
        raise ValueError(f"wd_thresh must be a non-negative float or int. Check water depth threshold argument ({_ARG_NAMES.WD_THRESH})")
        
    if bbox is not None:
        if len(bbox) != 4:
            raise ValueError(f"bbox must be a tuple of four floats (minx, miny, maxx, maxy). Check bbox argument ({_ARG_NAMES.BBOX})")
        if not all(isinstance(coord, (int, float)) for coord in bbox):
            raise TypeError(f"All coordinates in bbox must be int or float. Check bbox argument ({_ARG_NAMES.BBOX})")
        bbox = gpd.GeoDataFrame({'geometry': [box(*bbox)]}, crs="EPSG:4326")
    else:
        if buildings_filename is not None:
            buildings_gdf = gpd.read_file(buildings_filename)
            bbox = gpd.GeoDataFrame({'geometry': [box(*buildings_gdf.total_bounds)]}, crs=_utils.get_geodataframe_crs(buildings_gdf))
        else:
            bbox = gpd.GeoDataFrame({'geometry': [box(*_utils.get_raster_bounds(waterdepth_filename))]}, crs=_utils.get_raster_crs(waterdepth_filename))
    
    if out is not None:
        if type(out) is not str:
            raise TypeError(f"out must be a string. Check output file argument ({_ARG_NAMES.OUT})")
        _valid_out_format = {'geojson', 'json', 'shp', 'gpkg'}
        out_ext = filesystem.justext(out)
        if out_ext not in _valid_out_format:
            raise ValueError(f"out must be a file path ending with one of the following extensions: {', '.join(_valid_out_format)}. Check output file argument ({_ARG_NAMES.OUT})")
    else:
        out = os.path.join(os.getcwd(), f"safer_buildings_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.geojson")
        
    if t_srs is not None:
        if type(t_srs) is not str:
            raise TypeError(f"t_srs must be a string. Check target SRS argument ({_ARG_NAMES.T_SRS})")
        if not t_srs.startswith("EPSG:"):
            raise ValueError(f"t_srs must be an EPSG code (e.g., 'EPSG:4326'). Check target SRS argument ({_ARG_NAMES.T_SRS})")
        try:
            epsg_code = int(t_srs.split(":")[1])
        except ValueError:
            raise ValueError(f"Invalid EPSG code in t_srs. Check target SRS argument ({_ARG_NAMES.T_SRS})")
        srs = osr.SpatialReference()
        valid = srs.ImportFromEPSG(epsg_code)
        if valid != 0:
            raise ValueError(f"Invalid EPSG code: {epsg_code}. Check target SRS argument ({_ARG_NAMES.T_SRS})")
    else:
        if buildings_filename is not None:
            t_srs = _utils.get_geodataframe_crs(gpd.read_file(buildings_filename))
        else:
            t_srs = _utils.get_raster_crs(waterdepth_filename)
        
        
    if provider is None:
        if buildings_filename is None:
            provider = _consts._OVERTURE_PROVIDER
            Logger.debug(f"## No provider specified. Defaulting to {_consts._OVERTURE_PROVIDER}.")
    else:
        if type(provider) is not str:
            raise TypeError("provider must be a string")
        if provider.startswith(_consts._RER_REST_PROVIDER):
            Logger.debug(f"## Init {_consts._RER_REST_PROVIDER} before validating the provider")
            _consts.init_rer_rest_layers()
            if provider != _consts._RER_REST_PROVIDER:
                if len(provider.split('/')) < 2:
                    raise ValueError(f"{_consts._RER_REST_PROVIDER} provider must be '{_consts._RER_REST_PROVIDER}' or in the format '{_consts._RER_REST_PROVIDER}/<service_id>'. At least one service_id must be provided. MULTIPLE service_ids can be specified by '/' separated list, e.g. '{_consts._RER_REST_PROVIDER}/30/31/32'. Check provider argument ({_ARG_NAMES.PROVIDER})")
                service_ids = provider.split('/')[1:]
                for service_id in service_ids:
                    provider_service = f'{_consts._RER_REST_PROVIDER}/{service_id}'
                    if provider_service not in _consts._PROVIDERS:
                        raise ValueError(f"Invalid provider: {provider_service}. Valid providers are: {_consts._PROVIDERS}. Check provider argument ({_ARG_NAMES.PROVIDER})")
        elif provider.startswith(_consts._VENEZIA_WFS_PROVIDER):
            Logger.debug(f"## Init {_consts._VENEZIA_WFS_PROVIDER} before validating the provider")
            _consts.init_venezia_wfs_layers()
            if provider not in (_consts._VENEZIA_WFS_PROVIDER, _consts._VENEZIA_WFS_CRITICAL_SITES_PROVIDER):
                if len(provider.split('/')) < 2:
                    raise ValueError(f"{_consts._VENEZIA_WFS_PROVIDER} provider must be '{_consts._VENEZIA_WFS_PROVIDER}' or in the format '{_consts._VENEZIA_WFS_PROVIDER}/<feature_name>'. At least one feature_name must be provided. MULTIPLE feature_names can be specified by '/' separated list, e.g. '{_consts._VENEZIA_WFS_PROVIDER}/feature1/feature2'. Check provider argument ({_ARG_NAMES.PROVIDER})")
                feature_names = provider.split('/')[1:]
                for feature_name in feature_names:
                    provider_feature = f'{_consts._VENEZIA_WFS_PROVIDER}/{feature_name}'
                    if provider_feature not in _consts._PROVIDERS:
                        raise ValueError(f"Invalid provider: {provider_feature}. Valid providers are: {_consts._PROVIDERS}. Check provider argument ({_ARG_NAMES.PROVIDER})")
        elif provider not in _consts._PROVIDERS:
            raise ValueError(f"Invalid provider: {provider}. Valid providers are: {_consts._PROVIDERS}. Check provider argument ({_ARG_NAMES.PROVIDER})")
    
    if flood_mode is not None:
        if type(flood_mode) is not str:
            raise TypeError(f"flood_mode must be a string. Check flood_mode argument ({_ARG_NAMES.FLOOD_MODE})")
        flood_mode = flood_mode.upper()
        if flood_mode not in _consts._FLOOD_MODES:
            raise ValueError(f"Invalid flood mode: {flood_mode}. Valid flood modes are: {_consts._FLOOD_MODES}. Check flood_mode argument ({_ARG_NAMES.FLOOD_MODE})")
        if flood_mode == _consts.FloodModes.BUFFER and provider.startswith(_consts._VENEZIA_WFS_PROVIDER):
            raise ValueError(f'## Flood mode "{_consts.FloodModes.BUFFER}" could end in incorrect results as {_consts._VENEZIA_WFS_PROVIDER} provider provides also some open-areas. Consider using a {_ARG_NAMES.FLOOD_MODE} valued as "{_consts.FloodModes.IN_AREA}" or "{_consts.FloodModes.ALL}" instead.')
        if flood_mode in (_consts.FloodModes.IN_AREA, _consts._FLOOD_MODES.ALL) and (provider == _consts.OVERTURE or provider.startswith(_consts._RER_REST_PROVIDER)):
            raise ValueError(f'## Flood mode "{flood_mode}" could end in incorrect results as {_consts.OVERTURE} and {_consts._RER_REST_PROVIDER} providers provide only buildings. Consider using a {_ARG_NAMES.FLOOD_MODE} valued as "{_consts.FloodModes.BUFFER}" instead.')
    else:
        flood_mode = _consts.FloodModes.BUFFER
        if provider.startswith(_consts._VENEZIA_WFS_PROVIDER):
            flood_mode = _consts.FloodModes.ALL

    if feature_filters is not None:
        if not isinstance(feature_filters, (dict, list)):
            raise TypeError(f"filters must be a dictionary or a list. Check filters argument ({_ARG_NAMES.FILTERS})")
        if isinstance(feature_filters, dict):
            feature_filters = [feature_filters]
        for idxf, filters in enumerate(feature_filters):
            for f_key, f_value in filters.items():
                if not isinstance(f_key, str):
                    raise TypeError(f"Filter keys must be strings. Check filters argument ({_ARG_NAMES.FILTERS})")
                if not isinstance(f_value, (list, str, int)):
                    raise TypeError(f"Filter values must be a list, string, or integer. Check filters argument ({_ARG_NAMES.FILTERS})")
                if not isinstance(f_value, list):
                    f_value = [f_value]
                filters[f_key] = f_value
            feature_filters[idxf] = filters                            
    else:
        feature_filters = []
        
    if only_flood is None:
        only_flood = False
    if type(only_flood) is not bool:
        raise TypeError(f"only_flood must be a boolean value. Check only_flood argument ({_ARG_NAMES.ONLY_FLOOD})")
        
    if compute_stats is None:
        compute_stats = False
    if type(compute_stats) is not bool:
        raise TypeError(f"stats must be a boolean value. Check stats argument ({_ARG_NAMES.STATS})")
    
    if compute_summary is None:
        compute_summary = False
    if type(compute_summary) is not bool:
        raise TypeError(f"summary must be a boolean value. Check summary argument ({_ARG_NAMES.SUMMARY})")
    
    if summary_on is not None:
        if isinstance(summary_on, str):
            summary_on = [summary_on]
        if not isinstance(summary_on, list):
            raise TypeError(f"summary_on must be a string or a list of strings. Check summary_on argument ({_ARG_NAMES.SUMMARY_ON})")
        if not all([isinstance(item, str) for item in summary_on]):
            raise TypeError(f"All items in summary_on must be strings. Check summary_on argument ({_ARG_NAMES.SUMMARY_ON})")
        compute_summary = True

    if add_ops is not None:
        for op_name, op_args in add_ops.items():
            if op_name not in module_add_ops._ADD_OPS:
                raise ValueError(f"Invalid additional operation: {op_name}. Valid operations are: {list(module_add_ops._ADD_OPS.keys())}. Check add_ops argument ({_ARG_NAMES.ADD_OPS})")
            if not isinstance(op_args, dict):
                raise TypeError(f"Arguments for additional operation '{op_name}' must be a dictionary. Check add_ops argument ({_ARG_NAMES.ADD_OPS})")
            if provider is not None and provider.split('/')[0] not in module_add_ops._ADD_OPS[op_name]['providers']:
                raise ValueError(f"Additional operation '{op_name}' is not supported by provider '{provider}'. Supported providers are: {module_add_ops._ADD_OPS[op_name]['providers']}. Check add_ops argument ({_ARG_NAMES.ADD_OPS})")
            for argn,_ in op_args.items():
                if not isinstance(argn, str):
                    raise TypeError(f"Argument names for additional operation '{op_name}' must be strings. Check add_ops argument ({_ARG_NAMES.ADD_OPS})")
                if argn not in module_add_ops._ADD_OPS[op_name]['args']:
                    raise ValueError(f"Invalid argument '{argn}' for additional operation '{op_name}'. Valid arguments are: {module_add_ops._ADD_OPS[op_name]['args']}. Check add_ops argument ({_ARG_NAMES.ADD_OPS})")
        add_ops = [module_add_ops._ADD_OPS[op_name]['class'](**op_args) for op_name, op_args in add_ops.items()]
            
    
    if out_geojson is None:
        out_geojson = False
    if type(out_geojson) is not bool:
        raise TypeError(f"out_geojson must be a boolean value. Check out_geojson argument ({_ARG_NAMES.OUT_GEOJSON})")
    
        
    Logger.debug("## Input arguments validated successfully.")
    Logger.debug(f"### Water depth file: {waterdepth_filename}")
    Logger.debug(f"### Buildings file: {buildings_filename}")
    Logger.debug(f"### Water depth threshold: {wd_thresh}")
    Logger.debug(f"### Bounding box (total bounds): {bbox.total_bounds}")
    Logger.debug(f"### Output file: {out}")
    Logger.debug(f"### Target SRS: {t_srs}")
    Logger.debug(f"### Providers: {provider}")
    Logger.debug(f"### Flood mode: {flood_mode}")
    Logger.debug(f"### Feature filters: {feature_filters}")
    Logger.debug(f"### Only flood: {only_flood}")
    Logger.debug(f"### Compute stats: {compute_stats}")
    Logger.debug(f"### Compute summary: {compute_summary}")
    Logger.debug(f"### Summary on: {summary_on}")
    Logger.debug(f"### Additional operations: {[op.name for op in add_ops] if add_ops else None}")
    Logger.debug(f"### Output GeoJSON: {out_geojson}")
        
    return waterdepth_filename, buildings_filename, wd_thresh, bbox, out, t_srs, provider, flood_mode, feature_filters, only_flood, compute_stats, compute_summary, summary_on, add_ops, out_geojson

