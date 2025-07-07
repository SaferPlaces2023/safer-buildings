import geopandas as gpd
from shapely.geometry import MultiPolygon
from shapely import intersects

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
    
    Logger.debug(f"## Compute ring geometries around buildings (radius: {_consts._RING_BUFFER_M} meters) ...")
    radius_buffer = _consts._RING_BUFFER_M * (1 if _utils.crs_is_projected(f'EPSG:{buildings.crs.to_epsg()}') else 1e-5)
    buildings_rings = _utils.get_polygon_ring(buildings, radius_buffer)
    buildings['ring_geometry'] = buildings_rings.geometry
    
    Logger.debug("## Get buildings with intersection between ring geometries and water depth polygons ...")
    buildings['flood_bounds'] = buildings.ring_geometry.apply(lambda rg: MultiPolygon(polygons=waterdepth_mask.cx[rg.bounds[0]:rg.bounds[2], rg.bounds[1]:rg.bounds[3]].geometry.tolist()))
    buildings['is_flooded'] = buildings.apply(lambda b: intersects(b.ring_geometry, b.flood_bounds) if not b.flood_bounds.is_empty else False, axis=1)
    
    Logger.debug(f"## Found {len(buildings[buildings['is_flooded']])} flooded buildings out of {len(buildings)} total buildings.")
    
    return buildings