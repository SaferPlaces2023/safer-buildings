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


def compute_ring_geometry(
    buildings: gpd.GeoDataFrame,
    radius_buffer: float
) -> gpd.GeoDataFrame:
    """
    Compute ring geometries around the features in a GeoDataFrame.
    """
    Logger.debug(f"## Compute ring geometries around features (radius: {radius_buffer} meters) ...")
    radius_buffer = _consts._RING_BUFFER_M * (1 if _utils.crs_is_projected(f'EPSG:{buildings.crs.to_epsg()}') else 1e-5)
    buildings_rings = _utils.get_polygon_ring(buildings, radius_buffer)
    buildings['ring_geometry'] = buildings_rings.geometry
    return buildings


def compute_flood_area(
    waterdepth_mask: gpd.GeoDataFrame,
    buildings: gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:
    """
    Compute flood area by intersecting water depth polygons with buildings.
    """
    Logger.debug("## Get buildings with intersection between ring geometries and water depth polygons ...")
    
    wd_tree = STRtree(waterdepth_mask.geometry.values)

    def get_intersecting_multipolygon(ring_geom):
        candidates = wd_tree.geometries.take(wd_tree.query(ring_geom)).tolist()
        return MultiPolygon(polygons=[g for g in candidates if ring_geom.intersects(g)])
    
    buildings['flood_area'] = buildings['ring_geometry'].apply(get_intersecting_multipolygon)
    return buildings
    

def compute_flooded_buildings(
    buildings: gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:
    """
    Check if buildings are flooded by intersecting ring geometries with flood area.
    """
    Logger.debug("## Check if buildings are flooded by intersecting ring geometries with flood area ...")
    
    not_empty = ~buildings['flood_area'].is_empty
    intersects = buildings['ring_geometry'].intersects(buildings['flood_area'])
    buildings['is_flooded'] = not_empty & intersects

    Logger.debug(f"## Found {len(buildings[buildings['is_flooded']])} flooded buildings out of {len(buildings)} total buildings.")
    
    return buildings
        


def get_flooded_buildings(
    waterdepth_mask: gpd.GeoDataFrame | None,
    buildings: gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:
    
    """ 
    Get flooded buildings by intersecting water depth polygons with buildings.
    """
    
    buildings = _utils.ensure_geodataframe_crs(buildings, _utils.get_geodataframe_crs(waterdepth_mask))
    
    buildings = compute_ring_geometry(buildings, _consts._RING_BUFFER_M)
    
    buildings = compute_flood_area(waterdepth_mask, buildings)
    
    buildings = compute_flooded_buildings(buildings)

    return buildings