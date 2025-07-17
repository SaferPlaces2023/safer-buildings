import os
import datetime

from osgeo import osr

from shapely.geometry import box
import geopandas as gpd

from . import _utils
from . import _consts
from . import module_s3

from .module_log import Logger

class _ARG_NAMES():
    WATER = ['--water', '--waterdepth', '--wd', '--water_filename']
    BUILDINGS = ['--building', '--buildings', '--buildings_filename']
    WD_THRESH = ['--wd_thresh', '--thresh']
    BBOX = ['--bbox']
    OUT = ['--out']
    T_SRS = ['--t_srs']
    PROVIDER = ['--provider']
    FILTERS = ['--filters']
    ONLY_FLOOD = ['--only_flood']
    STATS = ['--stats']
    SUMMARY = ['--summary']
    SUMMARY_ON = ['--summary_on']
    OUT_GEOJSON = ['--out_geojson']


def validate_args(
    waterdepth_filename: str,
    buildings_filename: str | None = None,
    wd_thresh: float = 0.5,
    bbox: tuple[float, float, float, float] | None = None,
    out: str | None = None,
    t_srs: str | None = None,
    provider: list[str] | None = None,
    feature_filters: dict[str, dict] | None = None,
    only_flood: bool = False,
    compute_stats: bool = False,
    compute_summary: bool = False,
    summary_on: str | list[str] | None = None,
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
    
    if buildings_filename is not None:
        if type(buildings_filename) is not str:
            raise TypeError(f"buildings must be a string. Check buildings arguments ({_ARG_NAMES.BUILDINGS})")
        if not buildings_filename.startswith('s3://') and os.path.isfile(buildings_filename) is False:
            raise FileNotFoundError(f"Buildings file not found: {buildings_filename}. Check buildings arguments ({_ARG_NAMES.BUILDINGS})")
        if buildings_filename.startswith('s3://'):
            buildings_filename_tmp = _utils.temp_filename(ext = _utils.justext(buildings_filename), prefix='safer-buildings_buildings')
            if buildings_filename.endswith('.shp'):
                add_ext = ['.shx', '.dbf', '.prj', '.cpg']
                for ext in add_ext:
                    dl_status = module_s3.s3_download(uri = buildings_filename.replace('.shp', ext), fileout=buildings_filename_tmp.replace('.shp', ext))
                    if dl_status is None:
                        raise FileNotFoundError(f"Failed to download buildings files from S3: {buildings_filename}. Check buildings arguments ({_ARG_NAMES.BUILDINGS})")
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
        if not (out.endswith('.geojson') or out.endswith('.json')):
            raise ValueError(f"out must be a file path ending with .geojson or .json. Check output file argument ({_ARG_NAMES.OUT})")
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
            raise ValueError(f"A provider must be provided if buildings_filename is not given. Check provider argument ({_ARG_NAMES.PROVIDER}) or buildings argument ({_ARG_NAMES.BUILDINGS})")
    else:
        if type(provider) is not str:
            raise TypeError("provider must be a string")
        if provider.startswith('RER-REST'):
            if provider != 'RER-REST':
                if len(provider.split('/')) < 2:
                    raise ValueError(f"RER-REST provider must be 'RER-REST' or in the format 'RER-REST/<service_id>'. At least one service_id must be provided. MULTIPLE service_ids can be specified by '/' separated list, e.g. 'RER-REST/30/31/32'. Check provider argument ({_ARG_NAMES.PROVIDER})")
                service_ids = provider.split('/')[1:]
                for service_id in service_ids:
                    provider_service = f'RER-REST/{service_id}'
                    if provider_service not in _consts._PROVIDERS:
                        raise ValueError(f"Invalid provider: {provider_service}. Valid providers are: {_consts._PROVIDERS}. Check provider argument ({_ARG_NAMES.PROVIDER})")
        elif provider.startswith('VENEZIA-WFS'):
            if provider != 'VENEZIA-WFS':
                if len(provider.split('/')) < 2:
                    raise ValueError(f"VENEZIA-WFS provider must be 'VENEZIA-WFS' or in the format 'VENEZIA-WFS/<feature_name>'. At least one feature_name must be provided. MULTIPLE feature_names can be specified by '/' separated list, e.g. 'VENEZIA-WFS/feature1/feature2'. Check provider argument ({_ARG_NAMES.PROVIDER})")
                feature_names = provider.split('/')[1:]
                for feature_name in feature_names:
                    provider_feature = f'VENEZIA-WFS/{feature_name}'
                    if provider_feature not in _consts._PROVIDERS:
                        raise ValueError(f"Invalid provider: {provider_feature}. Valid providers are: {_consts._PROVIDERS}. Check provider argument ({_ARG_NAMES.PROVIDER})")
        elif provider not in _consts._PROVIDERS:
            raise ValueError(f"Invalid provider: {provider}. Valid providers are: {_consts._PROVIDERS}. Check provider argument ({_ARG_NAMES.PROVIDER})")
    
    
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
    Logger.debug(f"### Feature filters: {feature_filters}")
    Logger.debug(f"### Only flood: {only_flood}")
    Logger.debug(f"### Compute stats: {compute_stats}")
    Logger.debug(f"### Compute summary: {compute_summary}")
    Logger.debug(f"### Summary on: {summary_on}")
    Logger.debug(f"### Output GeoJSON: {out_geojson}")
        
    return waterdepth_filename, buildings_filename, wd_thresh, bbox, out, t_srs, provider, feature_filters, only_flood, compute_stats, compute_summary, summary_on, out_geojson

