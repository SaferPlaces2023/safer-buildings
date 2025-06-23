import geopandas as gpd

from . import _utils
from . import _consts

from .module_log import Logger


def get_flooded_buildings(
    waterdepth_mask: gpd.GeoDataFrame | None,
    buildings: gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:
    
    """ 
    Get flooded buildings by intersecting water depth polygons with buildings.
    """
    
    buildings = _utils.ensure_geodataframe_crs(buildings, _utils.get_geodataframe_crs(waterdepth_mask))
    
    radius_buffer = _consts._RING_BUFFER_M * (1 if _utils.crs_is_projected(f'EPSG:{buildings.crs.to_epsg()}') else 1e-5)
    
    # Intersect with water depth polygons
    buildings['__tmp_identifier__'] = buildings.index.to_list()
    buildings_rings = _utils.get_polygon_ring(buildings, radius_buffer)
    
    flooded_buildings = gpd.overlay(buildings_rings, waterdepth_mask, how='intersection')
    
    buildings['is_flooded'] = buildings['__tmp_identifier__'].apply(lambda tmp_id: tmp_id in flooded_buildings['__tmp_identifier__'].to_list())
    buildings.drop(columns=['__tmp_identifier__'], inplace=True)
    
    Logger.debug(f"### Found {len(flooded_buildings)} flooded buildings out of {len(buildings)} total buildings.")
    
    return buildings