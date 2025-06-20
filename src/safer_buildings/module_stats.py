import numpy as np
import pandas as pd
import geopandas as gpd

from osgeo import gdal

from . import _utils
from . import _consts

from .module_log import Logger


def filter_by_feature(
    gdf: gpd.GeoDataFrame,
    feature_filters: list[dict[str, list]],
    only_flood: bool
) -> gpd.GeoDataFrame:
    """
    Filter a GeoDataFrame by a list of filters.
    """

    if only_flood:
        Logger.info('## Filtering only flooded buildings ...')
        filtered_flooded_buildings = filtered_flooded_buildings[filtered_flooded_buildings['is_flooded']]
        Logger.info(f"### Only flooded buildings retained: {len(filtered_flooded_buildings)} out of {len(gdf)}.")

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
    

def compute_wd_stats(
    waterdepth_filename: str,
    waterdepth_mask: gpd.GeoDataFrame,
    buildings: gpd.GeoDataFrame
) -> gpd.GeoDataFrame:
    """
    Compute statistics on water depth around flooded buildings.
    """
    
    waterdepth_raster = gdal.Open(waterdepth_filename)

    def building_wd_stats(building_flood_area):
        flood_area_values = _utils.raster_sample_area(waterdepth_raster, building_flood_area)
        flood_area_stats = dict(pd.Series(flood_area_values).describe())
        return flood_area_stats 
        
    radius_buffer = _consts._RING_BUFFER_M * (1 if _utils.crs_is_projected(f'EPSG:{buildings.crs.to_epsg()}') else 1e-5)
    buildings_circles = buildings.buffer(radius_buffer)
    buildings_rings = buildings_circles.difference(buildings.geometry)
    builidngs_flood_area = buildings_rings.intersection(waterdepth_mask.geometry.iloc[0])

    flood_buildings_stats = [building_wd_stats(building_flood_area) if is_flooded else None for building_flood_area,is_flooded in zip(builidngs_flood_area, buildings['is_flooded'])]

    buildings['flood_wd_min'] = [stats['min'] if stats else None for stats in flood_buildings_stats]
    buildings['flood_wd_25perc'] = [stats['25%'] if stats else None for stats in flood_buildings_stats]
    buildings['flood_wd_mean'] = [stats['mean'] if stats else None for stats in flood_buildings_stats]
    buildings['flood_wd_median'] = [stats['50%'] if stats else None for stats in flood_buildings_stats]
    buildings['flood_wd_75perc'] = [stats['75%'] if stats else None for stats in flood_buildings_stats]
    buildings['flood_wd_max'] = [stats['max'] if stats else None for stats in flood_buildings_stats]

    return buildings


def compute_wd_summary(
    buildings: gpd.GeoDataFrame,
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
            'flooded_buildings': int(gdf['is_flooded'].sum()),
        }
        if include_stats:
            _base_summary.update({
                'flood_wd_min': float(np.nanmin(gdf['flood_wd_min'].values)),
                'flood_wd_25perc': float(np.nanpercentile(gdf['flood_wd_25perc'].values, 25)),
                'flood_wd_mean': float(np.nanmean(gdf['flood_wd_mean'].values)),
                'flood_wd_median': float(np.nanmedian(gdf['flood_wd_median'].values)),
                'flood_wd_75perc': float(np.nanpercentile(gdf['flood_wd_75perc'].values, 75)),
                'flood_wd_max': float(np.nanmax(gdf['flood_wd_max'].values))
            })
        _base_summary = {k: v if not np.isnan(v) else None for k, v in _base_summary.items()}
        return _base_summary
    
    summary['overall'] = base_summary(buildings)
    
    class_column = None
    if provider == 'OVERTURE':
        class_column = 'subtype'
    elif provider.startswith('RER-REST'):
        class_column = 'rer_class'
    else:
        raise ValueError(f"Provider '{provider}' is not supported for summary computation. Available providers are: {_consts.PROVIDERS}.")     # DOC: Should never happen, but just in case.
    buildings[class_column] = buildings[class_column].fillna('other')
    summary['classes'] = {
        class_name: base_summary(class_gdf)
        for class_name, class_gdf in buildings.groupby(class_column)
    }
    
    return summary
