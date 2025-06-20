import os
import datetime

from osgeo import osr

from shapely.geometry import box
import geopandas as gpd

from . import _utils
from . import _consts
from . import module_s3

from .module_log import Logger




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
    compute_summary: bool = False
) -> tuple[str, str | None, float, gpd.GeoDataFrame, str, str, str, list[dict[str, list]], bool, bool, bool]:
    
    """
    Validate the input arguments for the flooded buildings analysis.
    """
        
    if waterdepth_filename is None:
        raise ValueError("waterdepth_filename must be provided.")
    if type(waterdepth_filename) is not str:
        raise TypeError("waterdepth_filename must be a string.")
    if not waterdepth_filename.startswith('s3://') and os.path.isfile(waterdepth_filename) is False:
        raise FileNotFoundError(f"Water depth file not found: {waterdepth_filename}")
    if waterdepth_filename.startswith('s3://'):
        waterdepth_filename_tmp = _utils.justfname(waterdepth_filename)
        module_s3.s3_download(uri = waterdepth_filename, fileout = waterdepth_filename_tmp)
        waterdepth_filename = waterdepth_filename_tmp
    
    if buildings_filename is not None:
        if type(buildings_filename) is not str:
            raise TypeError("buildings_filename must be a string.")
        if not buildings_filename.startswith('s3://') and os.path.isfile(buildings_filename) is False:
            raise FileNotFoundError(f"Buildings file not found: {buildings_filename}")
        if buildings_filename.startswith('s3://'):
            buildings_filename_tmp = _utils.justfname(buildings_filename)
            module_s3.s3_download(uri = buildings_filename, fileout = buildings_filename_tmp)
            buildings_filename = buildings_filename_tmp
        
    if wd_thresh is None:
        wd_thresh = 0.5
    if type(wd_thresh) not in (int, float):
        raise TypeError("wd_thresh must be a float or int.")
    if wd_thresh < 0:
        raise ValueError("wd_thresh must be a non-negative float or int.")
        
    if bbox is not None:
        if len(bbox) != 4:
            raise ValueError("bbox must be a tuple of four floats (minx, miny, maxx, maxy).")
        if not all(isinstance(coord, (int, float)) for coord in bbox):
            raise TypeError("All coordinates in bbox must be int or float.")
        bbox = gpd.GeoDataFrame({'geometry': [box(*bbox)]}, crs="EPSG:4326")
    else:
        bbox = gpd.GeoDataFrame({'geometry': [box(*_utils.get_raster_bounds(waterdepth_filename))]}, crs=_utils.get_raster_crs(waterdepth_filename))
    
    if out is not None:
        if type(out) is not str:
            raise TypeError("out must be a string.") 
        if _utils.justext(out) not in ['.geojson', '.json']:
            raise ValueError("out must be a file path ending with .geojson or .json.")
    else:
        out = os.path.join(os.getcwd(), f"safer_buildings_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.geojson")
        
    if t_srs is not None:
        if type(t_srs) is not str:
            raise TypeError("t_srs must be a string.")
        if not t_srs.startswith("EPSG:"):
            raise ValueError("t_srs must be an EPSG code (e.g., 'EPSG:4326').")
        try:
            epsg_code = int(t_srs.split(":")[1])
        except ValueError:
            raise ValueError("Invalid EPSG code in t_srs.")
        srs = osr.SpatialReference()
        valid = srs.ImportFromEPSG(epsg_code)
        if valid != 0:
            raise ValueError(f"Invalid EPSG code: {epsg_code}")
    else:
        t_srs = _utils.get_raster_crs(waterdepth_filename)
        
        
    if provider is None:
        if buildings_filename is not None:
            raise ValueError("A provider must be provided if buildings_filename is given.")
        else:
            provider = 'OVERTURE'
    if type(provider) is not str:
        raise TypeError("provider must be a string")
    if provider.startswith('RER-REST'):
        if len(provider.split('/')) < 2:
            raise ValueError("RER-REST provider must be in the format 'RER-REST/<service_id>'. At least one service_id must be provided. MULTIPLE service_ids can be specified by '/' separated list, e.g. 'RER-REST/30/31/32'.")
        service_ids = provider.split('/')[1:]
        for service_id in service_ids:
            provider_service = f'RER-REST/{service_id}'
            if provider_service not in _consts._PROVIDERS:
                raise ValueError(f"Invalid provider: {provider_service}. Valid providers are: {_consts._PROVIDERS}.")
    elif provider not in _consts._PROVIDERS:
        raise ValueError(f"Invalid provider: {provider}. Valid providers are: {_consts._PROVIDERS}.")
    
    
    if feature_filters is not None:
        if not isinstance(feature_filters, (dict, list)):
            raise TypeError("feature_filters must be a dictionary or a list.")
        if isinstance(feature_filters, dict):
            feature_filters = [feature_filters]
        for idxf, filters in enumerate(feature_filters):
            for f_key, f_value in filters.items():
                if not isinstance(f_key, str):
                    raise TypeError("Filter keys must be strings.")
                if not isinstance(f_value, (list, str, int)):
                    raise TypeError("Filter values must be a list, string, or integer.")
                if not isinstance(f_value, list):
                    f_value = [f_value]
                filters[f_key] = f_value
            feature_filters[idxf] = filters                            
    else:
        feature_filters = []
        
    if only_flood is None:
        only_flood = False
    if type(only_flood) is not bool:
        raise TypeError("only_flood must be a boolean value.")
        
    if compute_stats is None:
        compute_stats = False
    if type(compute_stats) is not bool:
        raise TypeError("compute_stats must be a boolean value.")
    
    if compute_summary is None:
        compute_summary = False
    if type(compute_summary) is not bool:
        raise TypeError("compute_summary must be a boolean value.")
        
    Logger.info("## Input arguments validated successfully.")
    Logger.info(f"### Water depth file: {waterdepth_filename}")
    Logger.info(f"### Buildings file: {buildings_filename}")
    Logger.info(f"### Water depth threshold: {wd_thresh}")
    Logger.info(f"### Bounding box (total bounds): {bbox.total_bounds}")
    Logger.info(f"### Output file: {out}")
    Logger.info(f"### Target SRS: {t_srs}")
    Logger.info(f"### Providers: {provider}")
    Logger.info(f"### Feature filters: {feature_filters}")
    Logger.info(f"### Only flood: {only_flood}")
    Logger.info(f"### Compute stats: {compute_stats}")
    Logger.info(f"### Compute summary: {compute_summary}")
        
    return waterdepth_filename, buildings_filename, wd_thresh, bbox, out, t_srs, provider, feature_filters, only_flood, compute_stats, compute_summary

