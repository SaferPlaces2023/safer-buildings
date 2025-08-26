import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import MultiPolygon
from shapely import intersects, intersection

import xarray as xr
import rioxarray

from osgeo import gdal

from . import _utils
from . import _consts
from . import module_flood

from .module_log import Logger


def filter_by_feature(
    gdf: gpd.GeoDataFrame,
    feature_filters: list[dict[str, list]]
) -> gpd.GeoDataFrame:
    """
    Filter a GeoDataFrame by a list of filters.
    """

    if len(feature_filters) > 0:
        filtered_gdfs = []
        for filters in feature_filters:
            
            or_gdf = gdf.copy()
            
            for f_idx, (f_key, f_value) in enumerate(filters.items()):
                if f_key not in gdf.columns:
                    raise ValueError(f"Filter key '{f_key}' not found a valid feature. Avaliable features for this request are: {or_gdf.columns.tolist()}")
                and_gdf = or_gdf.copy() if f_idx == 0 else and_gdf
                and_gdf = and_gdf[and_gdf[f_key].isin(f_value)]
            filtered_gdfs.append(and_gdf)
        
        filtered_gdf = pd.concat(filtered_gdfs, ignore_index=True)
        filtered_gdf = filtered_gdf.drop_duplicates()
        
        return filtered_gdf
    else:
        return gdf
    

def filter_only_flood(
    gdf: gpd.GeoDataFrame,
    only_flood: bool
):
    """
    Filter a GeoDataFrame to only include flooded buildings.
    """
    
    if only_flood:
        Logger.debug('## Filtering only flooded buildings ...')
        gdf = gdf[gdf[_consts._COL_IS_FLOODED]].copy().reset_index(drop=True)
        Logger.debug(f"### Only flooded buildings retained: {len(gdf)} out of {len(gdf)}.")
    
    return gdf
    

def compute_wd_stats(
    waterdepth_filename: str,
    waterdepth_mask: gpd.GeoDataFrame,
    waterdepth_thresh: float,
    buildings: gpd.GeoDataFrame,
    flood_mode: str
) -> gpd.GeoDataFrame:
    """
    Compute statistics on water depth around flooded buildings.
    """

    # DOC: WD values as netcdf so we can access by coord indexing + filter on thresh
    waterdepth_ds = rioxarray.open_rasterio(waterdepth_filename).sel(band=1)
    waterdepth_ds = xr.where(waterdepth_ds > waterdepth_thresh, waterdepth_ds, np.nan).rio.write_crs(waterdepth_ds.rio.crs)
    
    if _consts._COL_FLOOD_ROI not in buildings.columns:
        buildings = module_flood.compute_flood_roi_geometry(buildings, flood_mode)
        
    if _consts._COL_FLOOD_AREA not in buildings.columns:
        buildings = module_flood.compute_flood_area(waterdepth_mask, buildings)
    
    if _consts._COL_IS_FLOODED not in buildings.columns:
        buildings = module_flood.compute_flooded_buildings(buildings)
    
    Logger.debug("## Compute intersection areas between building ring geometries and flood areas ...")    
    buildings[_consts._COL_FLOOD_GEOMETRY] = buildings[_consts._COL_FLOOD_ROI].intersection(buildings[_consts._COL_FLOOD_AREA])
    
    Logger.debug("## Extract sampling points from flood geometries ...")
    wd_res = np.array(waterdepth_ds.rio.resolution())
    wd_res = np.abs(wd_res / 1) 
    buildings[_consts._COL_FLOOD_COORDS] = buildings[_consts._COL_FLOOD_GEOMETRY].apply(lambda fg: _utils.coords_in_poly(fg, res=wd_res, poly_buffer=wd_res[0]) if not fg.is_empty else np.nan)    
        
    Logger.debug("## Compute water depth values at flood points + descriptive statistics ...")
    fpxs = buildings[_consts._COL_FLOOD_COORDS].apply(lambda pts: xr.DataArray(pts[:,0], dims=['loc']) if isinstance(pts, np.ndarray) else np.nan)
    fpys = buildings[_consts._COL_FLOOD_COORDS].apply(lambda pts: xr.DataArray(pts[:,1], dims=['loc']) if isinstance(pts, np.ndarray) else np.nan)
    flood_vals = [waterdepth_ds.sel(x=fpx, y=fpy, method='nearest').values if isinstance(fpx, xr.DataArray) else np.nan for fpx, fpy in zip(fpxs, fpys)]
    flood_stats = pd.Series(flood_vals).apply(lambda v: dict(pd.Series(v).describe()) if v is not np.nan else np.nan)
    
    buildings['flood_wd_min'] = flood_stats.apply(lambda s: s['min'] if isinstance(s, dict) else np.nan)
    buildings['flood_wd_25perc'] = flood_stats.apply(lambda s: s['25%'] if isinstance(s, dict) else np.nan)
    buildings['flood_wd_mean'] = flood_stats.apply(lambda s: s['mean'] if isinstance(s, dict) else np.nan)
    buildings['flood_wd_median'] = flood_stats.apply(lambda s: s['50%'] if isinstance(s, dict) else np.nan)
    buildings['flood_wd_75perc'] = flood_stats.apply(lambda s: s['75%'] if isinstance(s, dict) else np.nan)
    buildings['flood_wd_max'] = flood_stats.apply(lambda s: s['max'] if isinstance(s, dict) else np.nan)
    
    return buildings


def compute_wd_summary(
    buildings: gpd.GeoDataFrame,
    summary_on: list[str] | None,
    provider: str,
    include_stats: bool
) -> dict:
    """
    Compute summary statistics for flooded buildings.
    This function will return a dictionary with aggregated statistics based on building type and class.
    """
    
    summary = dict()
    
    def base_summary(gdf):
        _base_summary = {
            'total_buildings': len(gdf),
            'flooded_buildings': int(gdf[_consts._COL_IS_FLOODED].sum()),
        }
        if include_stats:
            _base_summary.update({
                'flood_wd_min': float(np.nanmin(gdf['flood_wd_min'].values)) if len(list(filter(pd.notna, gdf['flood_wd_min'].values))) > 0 else np.nan,
                'flood_wd_25perc': float(np.nanpercentile(gdf['flood_wd_25perc'].values, 25)) if len(list(filter(pd.notna, gdf['flood_wd_25perc'].values))) > 0 else np.nan,
                'flood_wd_mean': float(np.nanmean(gdf['flood_wd_mean'].values)) if len(list(filter(pd.notna, gdf['flood_wd_mean'].values))) > 0 else np.nan,
                'flood_wd_median': float(np.nanmedian(gdf['flood_wd_median'].values)) if len(list(filter(pd.notna, gdf['flood_wd_median'].values))) > 0 else np.nan,
                'flood_wd_75perc': float(np.nanpercentile(gdf['flood_wd_75perc'].values, 75)) if len(list(filter(pd.notna, gdf['flood_wd_75perc'].values))) > 0 else np.nan,
                'flood_wd_max': float(np.nanmax(gdf['flood_wd_max'].values)) if len(list(filter(pd.notna, gdf['flood_wd_max'].values))) > 0 else np.nan,
            })
        _base_summary = {k: v if not np.isnan(v) else None for k, v in _base_summary.items()}
        return _base_summary
    
    summary['overall'] = base_summary(buildings)
    
    class_column = None
    if summary_on is not None:
        summary_on = [attr for attr in summary_on if attr in buildings.columns]
        if len(summary_on) == 0:
            Logger.warning("## No valid attributes provided for summary computation. Summary will be computed for all flooded buildings.")
            return summary
        else:
            buildings['_summary_category'] = buildings[summary_on].apply(lambda row: [v for v in row.values if not pd.isna(v)], axis=1)
            buildings['_summary_category'] = buildings['_summary_category'].apply(lambda x: ' / '.join(x) if isinstance(x, list) and len(x) > 0 else np.nan)
            class_column = '_summary_category'
    elif provider is not None:
        if provider == _consts._OVERTURE_PROVIDER and 'subtype' in buildings.columns:
            class_column = 'subtype'
        elif provider.startswith(_consts._RER_REST_PROVIDER) and 'service_class' in buildings.columns:
            class_column = 'service_class'
        elif provider.startswith(_consts._VENEZIA_WFS_PROVIDER) and 'service_id' in buildings.columns:
            class_column = 'service_id'
        else:
            return summary
    else:
        return summary
    
    buildings[class_column] = buildings[class_column].fillna('other')
    summary['classes'] = {
        class_name: base_summary(class_gdf)
        for class_name, class_gdf in buildings.groupby(class_column)
    }
    
    return summary
