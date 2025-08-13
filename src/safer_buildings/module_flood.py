import pandas as pd
import geopandas as gpd
from shapely.geometry import MultiPolygon
from shapely import intersects
from shapely.strtree import STRtree

from . import _utils
from . import _consts

from .module_log import Logger


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
    )
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
        buildings = compute_ring_geometry(buildings, _consts._RING_BUFFER_M)

    # DOC: If flood_mode is IN_AREA, use the building geometries as flood ROI.
    elif flood_mode == _consts.FloodModes.IN_AREA:
        buildings[_consts._COL_FLOOD_ROI] = buildings.geometry

    # DOC: If flood_mode is ALL, use the buffered building geometries as flood ROI.
    elif flood_mode == _consts.FloodModes.ALL:
        radius_buffer = _consts._RING_BUFFER_M * (1 if _utils.crs_is_projected(f'EPSG:{buildings.crs.to_epsg()}') else 1e-5)
        buildings[_consts._COL_FLOOD_ROI] = buildings.geometry.buffer(radius_buffer)

    return buildings


def compute_ring_geometry(
    buildings: gpd.GeoDataFrame,
    radius_buffer: float
) -> gpd.GeoDataFrame:
    """
    Compute ring geometries around the features in a GeoDataFrame.
    """
    Logger.debug(f"## Compute ring geometries around features (radius: {radius_buffer} meters) ...")
    radius_buffer = radius_buffer * (1 if _utils.crs_is_projected(f'EPSG:{buildings.crs.to_epsg()}') else 1e-5)
    buildings_rings = _utils.get_polygon_ring(buildings, radius_buffer)
    buildings[_consts._COL_FLOOD_ROI] = buildings_rings.geometry
    return buildings


def compute_flood_area(
    waterdepth_gdf: gpd.GeoDataFrame,
    buildings: gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:
    """
    Compute flood area by intersecting water depth polygons with buildings.
    """
    Logger.debug("## Get buildings with intersection between ring geometries and water depth polygons ...")
    
    wd_tree = STRtree(waterdepth_gdf.to_crs(_consts._EPSG_UTMxx).geometry.values)

    building_wd_query = wd_tree.query(buildings[_consts._COL_FLOOD_ROI].to_crs(_consts._EPSG_UTMxx).values, predicate='intersects')
    buildings[_consts._COL_FLOOD_AREA] = gpd.GeoSeries([MultiPolygon(polygons=[]) for _ in range(len(buildings))], crs=_consts._EPSG_UTMxx)
    
    buildings_flood_area = pd.DataFrame(building_wd_query.T, columns=['bld_idxs', 'wd_idxs']).groupby('bld_idxs').agg(
        lambda wd_idxs: MultiPolygon(polygons=waterdepth_gdf.loc[wd_idxs].geometry.values)
    ).reset_index().rename(columns={'bld_idxs': 'bld_idx', 'wd_idxs': 'geometry'})
    buildings_flood_area = gpd.GeoDataFrame(buildings_flood_area, geometry='geometry', crs=_consts._EPSG_UTMxx)
    
    buildings.loc[buildings_flood_area['bld_idx'].to_list(), _consts._COL_FLOOD_AREA] = gpd.GeoSeries(
        index = buildings_flood_area['bld_idx'].to_list(),
        data = buildings_flood_area['geometry'].values,
    )

    buildings[_consts._COL_FLOOD_AREA] = buildings[_consts._COL_FLOOD_AREA].to_crs(buildings.crs)

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
    
    buildings = compute_flood_area(waterdepth_gdf, buildings)
    
    buildings = compute_flooded_buildings(buildings)

    return buildings