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
        Logger.debug('## Filtering only flooded buildings ...')
        filtered_flooded_buildings = filtered_flooded_buildings[filtered_flooded_buildings['is_flooded']]
        Logger.debug(f"### Only flooded buildings retained: {len(filtered_flooded_buildings)} out of {len(gdf)}.")

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
    waterdepth_thresh: float,
    buildings: gpd.GeoDataFrame
) -> gpd.GeoDataFrame:
    """
    Compute statistics on water depth around flooded buildings.
    """

    # DOC: WD values as netcdf so we can access by coord indexing + filter on thresh
    waterdepth_ds = rioxarray.open_rasterio(waterdepth_filename).sel(band=1)
    waterdepth_ds = xr.where(waterdepth_ds > waterdepth_thresh, waterdepth_ds, np.nan).rio.write_crs(waterdepth_ds.rio.crs)
    
    if 'ring_geometry' not in buildings.columns:
        Logger.debug(f"## Compute ring geometries around buildings (radius: {_consts._RING_BUFFER_M} meters) ...")
        radius_buffer = _consts._RING_BUFFER_M * (1 if _utils.crs_is_projected(f'EPSG:{buildings.crs.to_epsg()}') else 1e-5)
        buildings_rings = _utils.get_polygon_ring(buildings, radius_buffer)
        buildings['ring_geometry'] = buildings_rings.geometry
        
    if 'flood_bounds' not in buildings.columns:
        buildings['flood_bounds'] = buildings.ring_geometry.apply(lambda rg: MultiPolygon(polygons=waterdepth_mask.cx[rg.bounds[0]:rg.bounds[2], rg.bounds[1]:rg.bounds[3]].geometry.tolist()))
    
    if 'is_flooded' not in buildings.columns:
        Logger.debug("## Get buildings with intersection between ring geometries and water depth polygons ...")
        buildings['is_flooded'] = buildings.apply(lambda b: intersects(b.ring_geometry, b.flood_bounds) if not b.flood_bounds.is_empty else False, axis=1)
    
    Logger.debug("## Compute intersection areas between building ring geometries and flood areas ...")    
    buildings['flood_geometry'] = buildings.apply(lambda b: intersection(b.ring_geometry, b.flood_bounds) if b.is_flooded else b.flood_bounds, axis=1)
    
    Logger.debug("## Extract sampling points from flood geometries ...")
    wd_res = np.array(waterdepth_ds.rio.resolution())
    wd_res = np.abs(wd_res / 1) 
    buildings['flood_coords'] = buildings.flood_geometry.apply(lambda fg: _utils.coords_in_poly(fg, res=wd_res, poly_buffer=wd_res[0]) if not fg.is_empty else np.nan)    
        
    Logger.debug("## Compute water depth values at flood points + descriptive statistics ...")
    fpxs = buildings['flood_coords'].apply(lambda pts: xr.DataArray(pts[:,0], dims=['loc']) if isinstance(pts, np.ndarray) else np.nan)
    fpys = buildings['flood_coords'].apply(lambda pts: xr.DataArray(pts[:,1], dims=['loc']) if isinstance(pts, np.ndarray) else np.nan)
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
    if summary_on is not None:
        discarded_attrs = set(buildings.columns) - set(summary_on)
        if len(discarded_attrs) > 0:
            Logger.warning(f"## The following attributes will discarded from summary computation: {', '.join(discarded_attrs)}.")
        summary_on = [attr for attr in summary_on if attr in buildings.columns]
        if len(summary_on) == 0:
            Logger.warning("## No valid attributes provided for summary computation. Summary will be computed for all flooded buildings.")
            return summary
        else:
            buildings['_summary_category'] = buildings[summary_on].apply(lambda row: [v for v in row.values if not pd.isna(v)], axis=1)
            buildings['_summary_category'] = buildings['_summary_category'].apply(lambda x: ' / '.join(x) if isinstance(x, list) and len(x) > 0 else np.nan)
            class_column = '_summary_category'
    elif provider is not None:
        if provider == 'OVERTURE':
            class_column = 'subtype'
        elif provider.startswith('RER-REST'):
            class_column = 'service_class'
        elif provider.startswith('VENEZIA-WFS'):
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
