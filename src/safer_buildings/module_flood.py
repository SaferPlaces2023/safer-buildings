import numpy as np
import pandas as pd

from shapely.geometry import MultiPolygon, MultiPoint
from shapely import points
from shapely.strtree import STRtree
import geopandas as gpd

import rioxarray

from . import _utils
from . import _consts

from .module_log import Logger



def pointify_waterdepth(
    waterdepth_filename: str,
    wd_threshold: float,
    bbox: gpd.GeoDataFrame | None = None
):
    # DOC: Open wd tif as rioxarray dataset
    waterdepth_ds = rioxarray.open_rasterio(waterdepth_filename).sel(band=1)
    
    # DOC: Read valueas and corrds as numpy arrays
    waterdepth_vals = waterdepth_ds.data
    x_coords = waterdepth_ds.x.values
    y_coords = waterdepth_ds.y.values
    
    # DOC: Filter water depth values by threshold
    y_idxs, x_idxs = np.where(waterdepth_vals > wd_threshold)
    x_coords = x_coords[x_idxs]
    y_coords = y_coords[y_idxs]
    waterdepth_vals = waterdepth_vals[y_idxs, x_idxs]
    
    # DOC: Create a GeoDataFrame with water depth points
    gdf_wd = gpd.GeoDataFrame(
        {'wd_value': waterdepth_vals},
        geometry = points(x_coords, y_coords, waterdepth_vals),
        crs = _utils.get_raster_crs(waterdepth_filename)
    )
    
    # DOC: Clip GeoDataFrame to bounding box if provided
    if bbox is not None:
        bbox = _utils.ensure_geodataframe_crs(bbox, _utils.get_geodataframe_crs(gdf_wd)).total_bounds
        gdf_wd = gdf_wd.cx[bbox[0]:bbox[2], bbox[1]:bbox[3]]    
    
    return gdf_wd


def poligonyze_waterdepth(
    waterdepth_filename: str,
    wd_threshold: float,
    bbox: gpd.GeoDataFrame | None = None
) -> gpd.GeoDataFrame:
    waterdepth_polygonized = _utils.polygonize_raster_valid_data(
        raster_filename = waterdepth_filename,
        band = 1,
        mask_builder = lambda wd: wd > wd_threshold,    # Significant water depth threshold
        bbox = bbox
    ).reset_index(drop=True)
    return waterdepth_polygonized


def compute_flood_roi_geometry(
    buildings: gpd.GeoDataFrame,
    flood_mode: str
) -> gpd.GeoDataFrame:
    """
    Compute ring geometries around buildings based on flood mode.
    """

    # DOC: If flood_mode is BUFFER, use a default buffer radius to compute ring geometries as flood ROI.
    if flood_mode == _consts.FloodModes.BUFFER:
        buildings[_consts._COL_FLOOD_ROI] = _utils.get_polygon_ring(buildings, _consts._RING_BUFFER_M).geometry

    # DOC: If flood_mode is IN_AREA, use the building geometries as flood ROI.
    elif flood_mode == _consts.FloodModes.IN_AREA:
        buildings[_consts._COL_FLOOD_ROI] = buildings.geometry.to_crs(_consts._EPSG_UTMxx).buffer(0).to_crs(buildings.crs)

    # DOC: If flood_mode is ALL, use the buffered building geometries as flood ROI.
    elif flood_mode == _consts.FloodModes.ALL:
        buildings[_consts._COL_FLOOD_ROI] = buildings.geometry.to_crs(_consts._EPSG_UTMxx).buffer(_consts._RING_BUFFER_M).to_crs(buildings.crs)

    return buildings


def compute_poly_flood_area(
    waterdepth_gdf: gpd.GeoDataFrame,
    buildings: gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:
    """
    Compute flood area by intersecting water depth polygons with buildings.
    """
    Logger.debug("## Get buildings with intersection between ring geometries and water depth polygons ...")
    
    buildings[_consts._COL_FLOOD_AREA] = gpd.GeoSeries([MultiPolygon(polygons=[]) for _ in range(len(buildings))], crs=buildings.crs)

    wd_tree = STRtree(waterdepth_gdf.to_crs(crs=_consts._EPSG_UTMxx).geometry.values)
    building_wd_query = wd_tree.query(buildings[_consts._COL_FLOOD_ROI].to_crs(crs=_consts._EPSG_UTMxx).values, predicate='intersects')
    
    buildings_flood_area = pd.DataFrame(building_wd_query.T, columns=['bld_idxs', 'wd_idxs']).groupby('bld_idxs').agg(
        lambda wd_idxs: MultiPolygon(polygons=waterdepth_gdf.loc[wd_idxs].geometry.values)
    ).reset_index().rename(columns={'bld_idxs': 'bld_idx', 'wd_idxs': 'geometry'})
    buildings_flood_area = gpd.GeoDataFrame(buildings_flood_area, geometry='geometry', crs=waterdepth_gdf.crs).to_crs(crs=buildings.crs)
    
    buildings.loc[buildings_flood_area['bld_idx'].to_list(), _consts._COL_FLOOD_AREA] = gpd.GeoSeries(
        index = buildings_flood_area['bld_idx'].to_list(),
        data = buildings_flood_area['geometry'].values,
        crs = buildings.crs
    )

    return buildings


def compute_points_flood_area(
    waterdepth_gdf: gpd.GeoDataFrame,
    buildings: gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:
    """ c
    
    
    Compute flood area by intersecting water depth points with buildings. 
    """
    Logger.debug("## Get buildings with intersection between ring geometries and water depth points ...")
    
    buildings[_consts._COL_FLOOD_AREA] = gpd.GeoSeries([MultiPoint(points=[]) for _ in range(len(buildings))], crs=buildings.crs)
    
    wd_tree = STRtree(waterdepth_gdf.to_crs(crs=_consts._EPSG_UTMxx).geometry.values)
    building_wd_query = wd_tree.query(buildings[_consts._COL_FLOOD_ROI].to_crs(crs=_consts._EPSG_UTMxx).values, predicate='intersects')
    
    buildings_flood_area = pd.DataFrame(building_wd_query.T, columns=['bld_idxs', 'wd_idxs']).groupby('bld_idxs').agg(
        lambda wd_idxs: MultiPoint(points=waterdepth_gdf.loc[wd_idxs].geometry.values)
    ).reset_index().rename(columns={'bld_idxs': 'bld_idx', 'wd_idxs': 'geometry'})
    buildings_flood_area = gpd.GeoDataFrame(buildings_flood_area, geometry='geometry', crs=waterdepth_gdf.crs).to_crs(crs=buildings.crs)
    
    buildings.loc[buildings_flood_area['bld_idx'].to_list(), _consts._COL_FLOOD_AREA] = gpd.GeoSeries(
        index = buildings_flood_area['bld_idx'].to_list(),
        data = buildings_flood_area['geometry'].values,
        crs = buildings.crs
    )
    
    return buildings
    

def compute_flooded_buildings(
    buildings: gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:
    """
    Check if buildings are flooded by intersecting ring geometries with flood area.
    """
    Logger.debug("## Check if buildings are flooded by intersecting ring geometries with flood area ...")
    
    buildings[_consts._COL_IS_FLOODED] = ~buildings[_consts._COL_FLOOD_AREA].is_empty
    
    Logger.debug(f"## Found {len(buildings[buildings[_consts._COL_IS_FLOODED]])} flooded buildings out of {len(buildings)} total buildings.")
    
    return buildings
        


def get_flooded_buildings(
    waterdepth_gdf: gpd.GeoDataFrame | None,
    buildings: gpd.GeoDataFrame,
    flood_mode: str
) -> gpd.GeoDataFrame:
    
    """ 
    Get flooded buildings by intersecting water depth polygons with buildings.
    """
    
    buildings = _utils.ensure_geodataframe_crs(buildings, _utils.get_geodataframe_crs(waterdepth_gdf))
    
    buildings = compute_flood_roi_geometry(buildings, flood_mode)
    
    # buildings = compute_poly_flood_area(waterdepth_gdf, buildings)
    buildings = compute_points_flood_area(waterdepth_gdf, buildings)
    
    buildings = compute_flooded_buildings(buildings)

    return buildings