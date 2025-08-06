import os
import sys
import logging
import datetime
import tempfile

import numpy as np
import pandas as pd
from pandas.api.types import is_datetime64_any_dtype

from osgeo import gdal, osr, ogr

from shapely.wkt import loads
from shapely.ops import unary_union
from shapely.geometry import Point, MultiPoint, LineString, MultiLineString, box, Polygon, MultiPolygon
import geopandas as gpd

from ._consts import _GARBAGE_TEMP_FILES, collect_garbage_temp_file
from .module_log import Logger
from .module_version import get_version

_base_temp_dir = tempfile.gettempdir()
_module_temp_dir = os.path.join(_base_temp_dir, 'safer_buildings')
os.makedirs(_module_temp_dir, exist_ok=True)



class CustomException(Exception):
    """
    Custom exception class for safer_buildings module.
    """
    def __init__(self, message, status_code, *args):
        super().__init__(message, *args)
        self.message = message
        self.status_code = status_code

    def from_exception(cls, e, status_code=500):
        """
        Create a CustomException from an existing exception.
        """
        if isinstance(e, CustomException):
            return e
        else:
            return cls(str(e), status_code)



def process_cli_args(
    version: bool = False,
    debug: bool = False,
    verbose: bool = False,
):
    if version:
        Logger.setLevel(logging.INFO)  # Set to ERROR to avoid debug logs in version output
        Logger.info(f"safer-buildings v-{get_version()}")
        return {'version': get_version()}
        
    if debug:
        Logger.setLevel(logging.DEBUG)
    elif verbose:
        Logger.setLevel(logging.INFO)
    else:
        Logger.setLevel(logging.ERROR)



def temp_filename(ext, prefix='', add_to_garbage_collection=True):
    temp_filepath = os.path.join(_module_temp_dir, f'{prefix}__{datetime.datetime.now().strftime("%Y%m%d_%H%M%S.%f")}.{ext}')
    if add_to_garbage_collection:
        collect_garbage_temp_file(temp_filepath)
    return temp_filepath


def clean_temp_files(from_garbage_collection=True):
    """
    Cleans up temporary files collected for garbage collection.
    """
    if from_garbage_collection:
        n_files = len(_GARBAGE_TEMP_FILES)
        for fp in _GARBAGE_TEMP_FILES:
            try:
                if os.path.exists(fp):
                    os.remove(fp)
            except Exception as e:
                n_files -= 1
                Logger.error(f"### Error removing temporary file {fp}: {e}")
        _GARBAGE_TEMP_FILES.clear()
        Logger.debug(f"## Removed {n_files} temporary files from garbage collection.")
    else:
        tmp_fps = [os.path.join(_module_temp_dir, f) for f in os.listdir(_module_temp_dir) if os.path.isfile(os.path.join(_module_temp_dir, f))]
        n_files = len(tmp_fps)
        for fp in tmp_fps:
            try:
                if os.path.exists(fp):
                    os.remove(fp)
            except Exception as e:
                n_files -= 1
                Logger.error(f"### Error removing temporary file {fp}: {e}")
        Logger.debug(f"## Removed {n_files} temporary files from module temp directory: {_module_temp_dir}.")


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


def safe_json_df(df):
    df = df_dt_col_to_isoformat(df)
    df = df_nparr_col_to_value(df)
    return df

def df_dt_col_to_isoformat(df):
    for col in df.columns:
        if is_datetime64_any_dtype(df[col]):
            df[col] = df[col].apply(lambda x: x.isoformat() if pd.notnull(x) else None)
    return df

def df_nparr_col_to_value(df):
    for col in df.columns:
        if df[col].apply(lambda x: isinstance(x, np.ndarray) and x.ndim == 1 and x.size == 1).all():
            df[col] = df[col].apply(lambda x: x.item())
    return df



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
    

def buffer_points(gdf, buffer_meters):
    og_crs = gdf.crs
    gdf.to_crs(epsg=3857, inplace=True)
    gdf['geometry'] = gdf.geometry.apply(lambda g: g.buffer(buffer_meters) if type(g) in [Point, MultiPoint, LineString, MultiLineString] else g)
    gdf = gdf.to_crs(og_crs)
    return gdf


def get_polygon_ring(gdf: gpd.GeoDataFrame, ring_buffer):
    """
    Get the exterior ring of the first polygon in the GeoDataFrame.
    """
    gdf_rings = gdf.copy()
    gdf_rings.to_crs(epsg=3857, inplace=True)
    gdf_rings['geometry'] = gdf_rings.buffer(ring_buffer).geometry
    gdf_rings['geometry'] = gdf_rings['geometry'].difference(gdf.to_crs(epsg=3857).geometry)
    gdf_rings.to_crs(epsg=gdf.crs.to_epsg(), inplace=True)
    return gdf_rings


def coords_in_poly(poly, res, poly_buffer=None):
    if poly_buffer is not None:
        poly = poly.buffer(poly_buffer)
    res_x, res_y = (res, res) if isinstance(res, (int, float)) else res
    minx, miny, maxx, maxy = poly.bounds
    x_vals = np.arange(minx-res_x, maxx+res_x, res_x)
    y_vals = np.arange(miny-res_y, maxy+res_y, res_y)
    X, Y = np.meshgrid(x_vals, y_vals)
    cartesian_product = pd.Series(map(Point, np.column_stack([X.ravel(), Y.ravel()])))
    points = cartesian_product[poly.contains(cartesian_product)]
    coords = np.stack([[p.x,p.y] for p in points])
    return coords
    


def get_raster_crs(raster_filename):
    """
    Get EPSG string code of a raster.
    """
    dataset = gdal.Open(raster_filename, gdal.GA_ReadOnly)
    wkt = dataset.GetProjection()
    srs = osr.SpatialReference()
    srs.ImportFromWkt(wkt)
    srs.AutoIdentifyEPSG()
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
    gdf = gpd.GeoDataFrame({'geometry': list(gdf.geometry.iloc[0].geoms)}, crs = gdf.crs) # Explode MultiPolygon into Polygons
    
    if bbox is not None:
        bbox = ensure_geodataframe_crs(bbox, get_raster_crs(raster_filename)).total_bounds
        gdf = gdf.cx[bbox[0]:bbox[2], bbox[1]:bbox[3]]
    
    return gdf


def raster_sample_area(raster, area):
    if type(raster) is str:
        raster_ds = gdal.Open(raster)
    else:
        raster_ds = raster

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