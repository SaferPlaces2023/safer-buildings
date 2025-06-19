import os
import datetime
import tempfile

import numpy as np

from osgeo import gdal, osr, ogr

from shapely.wkt import loads
from shapely.geometry import box, Polygon
import geopandas as gpd

_base_temp_dir = tempfile.gettempdir()
_module_temp_dir = os.path.join(_base_temp_dir, 'safer_buildings')
os.makedirs(_module_temp_dir, exist_ok=True)


def temp_filename(ext, prefix=''):
    return os.path.join(_module_temp_dir, f'{prefix}__{datetime.datetime.now().strftime("%Y%m%d_%H%M%S.%f")}.{ext}')


def normpath(pathname):
    """
    normpath
    """
    if not pathname:
        return ""
    return os.path.normpath(pathname.replace("\\", "/")).replace("\\", "/")


def juststem(pathname):
    """
    juststem
    """
    pathname = os.path.basename(pathname)
    root, _ = os.path.splitext(pathname)
    return root

def justpath(pathname, n=1):
    """
    justpath
    """
    for _ in range(n):
        pathname, _ = os.path.split(normpath(pathname))
    if pathname == "":
        return "."
    return normpath(pathname)

def justfname(pathname):
    """
    justfname - returns the basename
    """
    return normpath(os.path.basename(normpath(pathname)))

def justext(pathname):
    """
    justext
    """
    pathname = os.path.basename(normpath(pathname))
    _, ext = os.path.splitext(pathname)
    return ext.lstrip(".")

def forceext(pathname, newext):
    """
    forceext
    """
    root, _ = os.path.splitext(normpath(pathname))
    pathname = root + ("." + newext if len(newext.strip()) > 0 else "")
    return normpath(pathname)



def eedem_bbox_2_shapely_bbox(eedem_bbox):
    """
    Converts eedem's bounding box to a shapely object. min_lat,max_lat,min_lon,max_lon to min_lon,min_lat,max_lon,max_lat
    """
    return box(eedem_bbox[2], eedem_bbox[0], eedem_bbox[3], eedem_bbox[1])

def shapely_bbox_2_eedem_bbox(shapely_bbox):
    """
    Converts shapely's bounding box to a eedem compliant list. min_lon,min_lat,max_lon,max_lat to min_lat,max_lat,min_lon,max_lon
    """
    bounds = shapely_bbox.bounds
    return [bounds[1], bounds[3], bounds[0], bounds[2]]


def crs_is_projected(epsg_string):
    """
    Check if the EPSG code is projected.
    """
    epsg_code = int(epsg_string.split(":")[-1])
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg_code)
    return srs.IsProjected()

    
def get_geodataframe_crs(geo_df):
    epsg_code = geo_df.crs.to_epsg()
    if epsg_code is None:
        raise ValueError("GeoDataFrame does not have a defined CRS.")
    return f"EPSG:{epsg_code}"

def ensure_geodataframe_crs(geo_df, epsg_string):
    """
    Ensures that the geodataframe has the specified CRS, if so converts it.
    """
    current_crs = geo_df.crs
    if current_crs == epsg_string:
        return geo_df
    else:
        epsg_code = int(epsg_string.split(":")[-1])
        srs = osr.SpatialReference()
        _ = srs.ImportFromEPSG(epsg_code)
        target_crs_wkt = srs.ExportToWkt()
        return geo_df.to_crs(crs=target_crs_wkt)
    

def get_polygon_ring(gdf: gpd.GeoDataFrame, ring_buffer):
    """
    Get the exterior ring of the first polygon in the GeoDataFrame.
    """
    gdf_rings = gdf.copy()
    gdf_rings['geometry'] = gdf_rings.buffer(ring_buffer).geometry
    gdf_rings = gpd.overlay(gdf_rings, gdf, how='difference').set_crs(gdf.crs)
    return gdf_rings
    


def get_raster_crs(raster_filename):
    """
    Get EPSG string code of a raster.
    """
    dataset = gdal.Open(raster_filename, gdal.GA_ReadOnly)
    wkt = dataset.GetProjection()
    srs = osr.SpatialReference()
    srs.ImportFromWkt(wkt)
    epsg_string = f"EPSG:{srs.GetAuthorityCode(None)}"
    dataset = None 
    return epsg_string


def get_raster_bounds(raster_filename):
    """
    Get raster bounds in the form (xmin, ymin, xmax, ymax)
    """
    dataset = gdal.Open(raster_filename, gdal.GA_ReadOnly)
    
    # Ottiene la trasformazione affine
    transform = dataset.GetGeoTransform()
    width = dataset.RasterXSize
    height = dataset.RasterYSize

    # Calcola i limiti spaziali
    minx = transform[0]
    maxx = minx + transform[1] * width
    maxy = transform[3]
    miny = maxy + transform[5] * height

    dataset = None  # Chiude il file raster

    return minx, miny, maxx, maxy


def polygonize_raster_valid_data(raster_filename, band=1, mask_builder=None, bbox=None):
    dataset = gdal.Open(raster_filename)
    band = dataset.GetRasterBand(1)
    nodata_value = band.GetNoDataValue()
    raster_array = band.ReadAsArray()
    if mask_builder is None:
        mask = np.where(np.isnan(raster_array) | (raster_array == nodata_value), 0, 1)
    else:
        mask = np.where(np.isnan(raster_array) | (raster_array == nodata_value) | (~ mask_builder(raster_array)), 0, 1)

    mem_driver = gdal.GetDriverByName("MEM")
    mask_ds = mem_driver.Create("", dataset.RasterXSize, dataset.RasterYSize, 1, gdal.GDT_Byte)
    mask_ds.SetGeoTransform(dataset.GetGeoTransform())
    mask_ds.SetProjection(dataset.GetProjection())
    mask_ds.GetRasterBand(1).WriteArray(mask)

    vector_driver = ogr.GetDriverByName("Memory")
    vector_ds = vector_driver.CreateDataSource("")
    layer = vector_ds.CreateLayer("mask_layer", None, ogr.wkbPolygon)
    field_defn = ogr.FieldDefn("value", ogr.OFTInteger)
    layer.CreateField(field_defn)

    gdal.Polygonize(mask_ds.GetRasterBand(1), None, layer, 0, [], None)

    valid_data_polygons = []
    for feature in layer:
        if feature.GetField("value") == 1:
            geom = feature.GetGeometryRef()
            if geom is not None and geom.ExportToWkt():
                valid_data_polygons.append(loads(geom.ExportToWkt()))
    dataset = None
    # Unisci i poligoni in una singola geometria
    valid_area = gpd.GeoSeries(valid_data_polygons).union_all()
    gdf = gpd.GeoDataFrame({'geometry':[valid_area]}, crs=get_raster_crs(raster_filename))
    
    if bbox is not None:
        bbox = ensure_geodataframe_crs(bbox, get_raster_crs(raster_filename))
        gdf = gpd.overlay(gdf, bbox, how='intersection')
    
    return gdf


def raster_sample_area(raster_filename, area):
    raster_ds = gdal.Open(raster_filename)
    raster_proj = raster_ds.GetProjection()
    
    gt = raster_ds.GetGeoTransform()
    cols = raster_ds.RasterXSize
    rows = raster_ds.RasterYSize

    # Create in-memory raster to rasterize the polygon
    mem_driver = gdal.GetDriverByName('MEM')
    mask_ds = mem_driver.Create('', cols, rows, 1, gdal.GDT_Byte)
    mask_ds.SetGeoTransform(gt)
    mask_ds.SetProjection(raster_proj)

    # Convert GeoPandas geometry to OGR geometry
    ogr_polygon = ogr.CreateGeometryFromWkt(area.wkt)

    # Create layer in memory and rasterize
    mem_layer = ogr.GetDriverByName("Memory").CreateDataSource("")
    layer = mem_layer.CreateLayer("poly", None, ogr.wkbPolygon)
    feature = ogr.Feature(layer.GetLayerDefn())
    feature.SetGeometry(ogr_polygon)
    layer.CreateFeature(feature)

    # Rasterize
    gdal.RasterizeLayer(mask_ds, [1], layer, burn_values=[1])
    mask_band = mask_ds.GetRasterBand(1)
    mask_array = mask_band.ReadAsArray()

    # Read raster band (e.g., band 1)
    raster_band = raster_ds.GetRasterBand(1)
    raster_array = raster_band.ReadAsArray()

    # Apply the mask
    masked_values = raster_array[mask_array == 1]

    # Remove no-data values if necessary
    nodata = raster_band.GetNoDataValue()
    if nodata is not None:
        masked_values = masked_values[masked_values != nodata]
        
    return masked_values