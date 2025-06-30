import os
import sys
import json
import click
import time
import logging
import requests
import datetime

import numpy as np
import pandas as pd

from concurrent.futures import ProcessPoolExecutor

from osgeo import gdal, ogr, osr
from shapely.wkt import loads
from shapely.geometry import box, Point, Polygon, MultiPolygon, LineString, MultiLineString
import geopandas as gpd

import leafmap

from . import _utils
from .module_log import Logger
from . import module_logo, module_args, module_retriever, module_flood, module_stats, module_s3


from dotenv import load_dotenv
load_dotenv()  

# logging.basicConfig(format="%(asctime)s [%(levelname)s] %(message)s", stream=sys.stderr)
# logger = logging.getLogger(__name__)    
# logger.setLevel(logging.INFO)



# TODO: --list-providers → shows detailed description of providers
# TODO: --summary → add metadata with aggregate stats ( n buildings, n flooded buildings, min,avg,max wd for each category (provider dependend) of buildings)





# DOC: Aux functions


# def validate_args(
#     waterdepth_filename: str,
#     buildings_filename: str | None = None,
#     wd_thresh: float = 0.5,
#     bbox: tuple[float, float, float, float] | None = None,
#     out: str | None = None,
#     t_srs: str | None = None,
#     provider: list[str] | None = None,
#     feature_filters: dict[str, dict] | None = None,
#     only_flood: bool = False,
#     compute_stats: bool = False,
#     compute_summary: bool = False
# ) -> tuple[str, str | None, float, gpd.GeoDataFrame, str, str, str, list[dict[str, list]], bool, bool, bool]:
    
#     """
#     Validate the input arguments for the flooded buildings analysis.
#     """
        
#     if waterdepth_filename is None:
#         raise ValueError("waterdepth_filename must be provided.")
#     if type(waterdepth_filename) is not str:
#         raise TypeError("waterdepth_filename must be a string.")
#     if os.path.isfile(waterdepth_filename) is False:
#         raise FileNotFoundError(f"Water depth file not found: {waterdepth_filename}")
    
#     if buildings_filename is not None:
#         if type(buildings_filename) is not str:
#             raise TypeError("buildings_filename must be a string.")
#         if not buildings_filename.startswith('s3://') and os.path.isfile(buildings_filename) is False:
#             raise FileNotFoundError(f"Buildings file not found: {buildings_filename}")
        
#     if wd_thresh is None:
#         wd_thresh = 0.5
#     if type(wd_thresh) not in (int, float):
#         raise TypeError("wd_thresh must be a float or int.")
#     if wd_thresh < 0:
#         raise ValueError("wd_thresh must be a non-negative float or int.")
        
#     if bbox is not None:
#         if len(bbox) != 4:
#             raise ValueError("bbox must be a tuple of four floats (minx, miny, maxx, maxy).")
#         if not all(isinstance(coord, (int, float)) for coord in bbox):
#             raise TypeError("All coordinates in bbox must be int or float.")
#         bbox = gpd.GeoDataFrame({'geometry': [box(*bbox)]}, crs="EPSG:4326")
#     else:
#         bbox = gpd.GeoDataFrame({'geometry': [box(*utils.get_raster_bounds(waterdepth_filename))]}, crs=utils.get_raster_crs(waterdepth_filename))
    
#     if out is not None:
#         if type(out) is not str:
#             raise TypeError("out must be a string.") 
#     else:
#         out = os.path.join(os.getcwd(), f"safer_buildings_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.geojson")
        
#     if t_srs is not None:
#         if type(t_srs) is not str:
#             raise TypeError("t_srs must be a string.")
#         if not t_srs.startswith("EPSG:"):
#             raise ValueError("t_srs must be an EPSG code (e.g., 'EPSG:4326').")
#         try:
#             epsg_code = int(t_srs.split(":")[1])
#         except ValueError:
#             raise ValueError("Invalid EPSG code in t_srs.")
#         srs = osr.SpatialReference()
#         valid = srs.ImportFromEPSG(epsg_code)
#         if valid != 0:
#             raise ValueError(f"Invalid EPSG code: {epsg_code}")
#     else:
#         t_srs = utils.get_raster_crs(waterdepth_filename)
        
        
#     if provider is None:
#         if buildings_filename is not None:
#             raise ValueError("A provider must be provided if buildings_filename is given.")
#         else:
#             provider = 'OVERTURE'
#     if type(provider) is not str:
#         raise TypeError("provider must be a string")
#     if provider.startswith('RER-REST'):
#         if len(provider.split('/')) < 2:
#             raise ValueError("RER-REST provider must be in the format 'RER-REST/<service_id>'. At least one service_id must be provided. MULTIPLE service_ids can be specified by '/' separated list, e.g. 'RER-REST/30/31/32'.")
#         service_ids = provider.split('/')[1:]
#         for service_id in service_ids:
#             provider_service = f'RER-REST/{service_id}'
#             if provider_service not in _PROVIDERS:
#                 raise ValueError(f"Invalid provider: {provider_service}. Valid providers are: {_PROVIDERS}.")
#     elif provider not in _PROVIDERS:
#         raise ValueError(f"Invalid provider: {provider}. Valid providers are: {_PROVIDERS}.")
    
    
#     if feature_filters is not None:
#         if not isinstance(feature_filters, (dict, list)):
#             raise TypeError("feature_filters must be a dictionary or a list.")
#         if isinstance(feature_filters, dict):
#             feature_filters = [feature_filters]
#         for idxf, filters in enumerate(feature_filters):
#             for f_key, f_value in filters.items():
#                 if not isinstance(f_key, str):
#                     raise TypeError("Filter keys must be strings.")
#                 if not isinstance(f_value, (list, str, int)):
#                     raise TypeError("Filter values must be a list, string, or integer.")
#                 if not isinstance(f_value, list):
#                     f_value = [f_value]
#                 filters[f_key] = f_value
#             feature_filters[idxf] = filters                            
#     else:
#         feature_filters = []
        
#     if only_flood is None:
#         only_flood = False
#     if type(only_flood) is not bool:
#         raise TypeError("only_flood must be a boolean value.")
        
#     if compute_stats is None:
#         compute_stats = False
#     if type(compute_stats) is not bool:
#         raise TypeError("compute_stats must be a boolean value.")
    
#     if compute_summary is None:
#         compute_summary = False
#     if type(compute_summary) is not bool:
#         raise TypeError("compute_summary must be a boolean value.")
#     # if compute_summary:
#     #     compute_stats = True  # If summary is requested, stats must be computed as well.
        
#     Logger.info("## Input arguments validated successfully.")
#     Logger.info(f"### Water depth file: {waterdepth_filename}")
#     Logger.info(f"### Buildings file: {buildings_filename}")
#     Logger.info(f"### Water depth threshold: {wd_thresh}")
#     Logger.info(f"### Bounding box (total bounds): {bbox.total_bounds}")
#     Logger.info(f"### Output file: {out}")
#     Logger.info(f"### Target SRS: {t_srs}")
#     Logger.info(f"### Providers: {provider}")
#     Logger.info(f"### Feature filters: {feature_filters}")
#     Logger.info(f"### Only flood: {only_flood}")
#     Logger.info(f"### Compute stats: {compute_stats}")
#     Logger.info(f"### Compute summary: {compute_summary}")
        
#     return waterdepth_filename, buildings_filename, wd_thresh, bbox, out, t_srs, provider, feature_filters, only_flood, compute_stats, compute_summary


# def retrieve_buildings(
#     buildings_filename: str | None, 
#     bbox: tuple[float, float, float, float],
#     provider: str
# ) -> gpd.GeoDataFrame:
    
#     """
#     Retrieve buildings data from specified providers.
#     If buildings_filename is provided, it will be used directly.
#     If not, it will download buildings data from the specified providers.
#     """
    
#     def retrieve_overture():
#         columns = ["id", "geometry", "height", "subtype", "class", "is_underground"]
#         provider_buildings = leafmap.get_overture_data(
#             overture_type = "building", 
#             bbox = bbox.to_crs(epsg=4326).total_bounds.tolist(), 
#             columns = columns
#         )
#         return provider_buildings
    
    
#     if buildings_filename is not None:
#         provider_buildings = gpd.read_file(buildings_filename)
#         bbox = bbox.to_crs(utils.get_geodataframe_crs(provider_buildings))
#         provider_buildings = gpd.overlay(provider_buildings, bbox, how='intersection')
#         Logger.info(f"## Using provided buildings data from {buildings_filename}. Found {len(provider_buildings)} buildings.")
    
#     else:
#         Logger.info(f"## Retrieving buildings data from provider: {provider} ...")
            
#         buildings_filename = utils.temp_filename(ext='shp', prefix=f"{provider.replace('/','-')}_buildings")
        
#         if provider == 'OVERTURE':
#             provider_buildings = retrieve_overture()
            
#         elif provider.startswith('RER-REST'):
#             # DOC: Retrieve from all subservices of the requested one.
#             service_ids = list(map(int, provider.split('/')[1:]))
#             while any([RegioneEmiliaRomagnaLayers[RegioneEmiliaRomagnaLayers.id == service_id].iloc[0].subLayerIds is not None for service_id in service_ids]):
#                 for service_id in service_ids:
#                     sub_layers = RegioneEmiliaRomagnaLayers[RegioneEmiliaRomagnaLayers.id == service_id].iloc[0].subLayerIds
#                     if sub_layers is not None:
#                         service_ids.remove(service_id)
#                         service_ids.extend(sub_layers)
#                         service_ids = list(set(service_ids))
                        
#             def rest_service_retrieve(service_id):
#                 url = f"https://servizigis.regione.emilia-romagna.it/geoags/rest/services/portale/saferplaces/MapServer/{service_id}/query"
#                 bbox4326 = bbox.to_crs("EPSG:4326").geometry.total_bounds
#                 params = {
#                     "geometry": ','.join([str(b) for b in bbox4326]),
#                     "geometryType": "esriGeometryEnvelope",
#                     "inSR": "4326",
#                     "outSR": "4326",
#                     "spatialRel": "esriSpatialRelIntersects",
#                     "where": "1=1",
#                     "outFields": "*",
#                     "f": "json",
#                     "returnGeometry": "true"
#                 }
#                 headers = {
#                     "User-Agent": "QGIS",
#                     "Referer": "http://localhost"
#                 }
#                 response = requests.get(url, params=params, headers=headers)
#                 data = response.json()
                
#                 features =[f['attributes'] for f in data['features']]
                
#                 geometry_type = RegioneEmiliaRomagnaLayers[RegioneEmiliaRomagnaLayers.id == service_id].iloc[0].geometryType
#                 if geometry_type == 'esriGeometryPoint':
#                     geometries = [
#                         Point(f['geometry']['x'], f['geometry']['y']) 
#                         for f in data['features']
#                     ]
#                 elif geometry_type == 'esriGeometryPolygon':
#                     geometries = [
#                         MultiPolygon(
#                             [Polygon([point for point in ring ])
#                             for ring in f['geometry']['rings']]
#                         ) 
#                         for f in data['features']
#                     ]
#                 elif geometry_type == 'esriGeometryPolyline':
#                     geometries = [
#                         MultiLineString(
#                             [LineString([point for point in path ])
#                             for path in f['geometry']['paths']]
#                         ) 
#                         for f in data['features']
#                     ]
#                 else:
#                     raise ValueError(f"Unsupported geometry type for service {service_id}: {RegioneEmiliaRomagnaLayers[RegioneEmiliaRomagnaLayers.id == service_id].iloc[0].geometryType}")

#                 rest_gdf = gpd.GeoDataFrame(features, geometry=geometries, crs=f"EPSG:{data['spatialReference']['wkid']}")
#                 rest_gdf['rer_id'] = service_id
#                 rest_gdf['rer_class'] = RegioneEmiliaRomagnaLayers[RegioneEmiliaRomagnaLayers.id == service_id].name.iloc[0]
                    
#                 return rest_gdf
            
#             def overture_intersection(gdf_re):
#                 Logger.info(f"### Overture intersection for {len(gdf_re)} RER-REST Points.")
#                 gdf_ot = retrieve_overture()
#                 gdf_re['ot_id'] = gdf_re.geometry.apply(lambda geom: gdf_ot[gdf_ot.geometry.contains(geom)].id.values.tolist() if type(geom) is Point else None)
#                 gdf_re['ot_id'] = gdf_re['ot_id'].apply(lambda ids: ids[0] if ids is not None and len(ids) > 0 else None)
#                 gdf_re['geometry'] = gdf_re.apply(lambda row: gdf_ot[gdf_ot.id == row.ot_id].iloc[0].geometry if row.ot_id is not None else row.geometry, axis=1)
#                 Logger.info(f'### Overture intersection: Taking overture building from {len(gdf_re[gdf_re.ot_id.notnull()])} original RER-REST Points.')
#                 return gdf_re
    
#             def buffer_points(gdf, buffer_meters = 20):
#                 gdf = gdf.to_crs("EPSG:7791")
#                 gdf['geometry'] = gdf.geometry.apply(lambda g: g.buffer(buffer_meters) if type(g) in [Point, LineString, MultiLineString] else g)
#                 gdf = gdf.to_crs("EPSG:4326")
#                 return gdf
            
#             provider_buildings = pd.concat([rest_service_retrieve(service_id) for service_id in service_ids], ignore_index=True)
#             provider_buildings = overture_intersection(provider_buildings)
#             provider_buildings = buffer_points(provider_buildings)
            
#             provider_buildings.to_file(buildings_filename, driver='ESRI Shapefile', index=False)
        
#         else:
#             raise ValueError(f"Provider '{provider}' is not supported. Available providers are: {_PROVIDERS}.")

        
#         Logger.info(f"### Buildings data retrieved from {provider} saved at {buildings_filename}. (Found {len(provider_buildings)} buildings)")        
        
#     return provider_buildings


# def get_waterdepth_mask(
#     waterdepth_filename: str,
#     mask_builder: callable = lambda wd: wd > 0.5,  # Default mask for significant water depth
# ) -> gpd.GeoDataFrame:
    
#     """
#     Create a mask from the water depth raster file.
#     The mask is created by polygonizing the raster data based on the provided mask_builder function.
#     """
    
#     waterdepth_polygonized = _utils.polygonize_raster_valid_data(
#         raster_filename=waterdepth_filename,
#         band=1,
#         mask_builder=mask_builder
#     )
    
#     return waterdepth_polygonized


# def get_flooded_buildings(
#     waterdepth_mask: gpd.GeoDataFrame | None,
#     buildings: gpd.GeoDataFrame,
# ) -> gpd.GeoDataFrame:
    
#     """ 
#     Get flooded buildings by intersecting water depth polygons with buildings.
#     """
    
#     buildings = _utils.ensure_geodataframe_crs(buildings, _utils.get_geodataframe_crs(waterdepth_mask))
    
#     radius_buffer = _RING_BUFFER_M * (1 if _utils.crs_is_projected(f'EPSG:{buildings.crs.to_epsg()}') else 1e-5)
    
#     # Intersect with water depth polygons
#     buildings['__tmp_identifier__'] = buildings.index.to_list()
#     buildings_rings = _utils.get_polygon_ring(buildings, radius_buffer)
    
#     flooded_buildings = gpd.overlay(buildings_rings, waterdepth_mask, how='intersection')
    
#     buildings['is_flooded'] = buildings['__tmp_identifier__'].apply(lambda tmp_id: tmp_id in flooded_buildings['__tmp_identifier__'].to_list())
#     buildings.drop(columns=['__tmp_identifier__'], inplace=True)
    
#     Logger.info(f"### Found {len(flooded_buildings)} flooded buildings out of {len(buildings)} total buildings.")
    
#     return buildings


# def filter_by_feature(
#     gdf: gpd.GeoDataFrame,
#     feature_filters: list[dict[str, list]]
# ) -> gpd.GeoDataFrame:
#     """
#     Filter a GeoDataFrame by a list of filters.
#     """
#     if len(feature_filters) > 0:
#         filtered_gdfs = []
#         for filters in feature_filters:
            
#             or_gdf = gdf.copy()
            
#             for f_idx, (f_key, f_value) in enumerate(filters.items()):
#                 if f_key not in gdf.columns:
#                     raise ValueError(f"Filter key '{f_key}' not found a valid feature. Avaliable features for this request are: {or_gdf.columns.tolist()}")
#                 and_gdf = or_gdf.copy() if f_idx == 0 else and_gdf
#                 and_gdf = and_gdf[and_gdf[f_key].isin(f_value)]
#             filtered_gdfs.append(and_gdf)
        
#         filtered_gdf = pd.concat(filtered_gdfs, ignore_index=True)
#         filtered_gdf = filtered_gdf.drop_duplicates()
        
#         return filtered_gdf
#     else:
#         return gdf
    

# def compute_wd_stats(
#     waterdepth_filename: str,
#     waterdepth_mask: gpd.GeoDataFrame,
#     buildings: gpd.GeoDataFrame
# ) -> gpd.GeoDataFrame:
#     """
#     Compute statistics on water depth around flooded buildings.
#     """
    
#     waterdepth_raster = gdal.Open(waterdepth_filename)

#     def building_wd_stats(building_flood_area):
#         flood_area_values = _utils.raster_sample_area(waterdepth_raster, building_flood_area)
#         flood_area_stats = dict(pd.Series(flood_area_values).describe())
#         return flood_area_stats 
        
#     radius_buffer = _RING_BUFFER_M * (1 if _utils.crs_is_projected(f'EPSG:{buildings.crs.to_epsg()}') else 1e-5)
#     buildings_circles = buildings.buffer(radius_buffer)
#     buildings_rings = buildings_circles.difference(buildings.geometry)
#     builidngs_flood_area = buildings_rings.intersection(waterdepth_mask.geometry.iloc[0])

#     flood_buildings_stats = [building_wd_stats(building_flood_area) if is_flooded else None for building_flood_area,is_flooded in zip(builidngs_flood_area, buildings['is_flooded'])]

#     buildings['flood_wd_min'] = [stats['min'] if stats else None for stats in flood_buildings_stats]
#     buildings['flood_wd_25perc'] = [stats['25%'] if stats else None for stats in flood_buildings_stats]
#     buildings['flood_wd_mean'] = [stats['mean'] if stats else None for stats in flood_buildings_stats]
#     buildings['flood_wd_median'] = [stats['50%'] if stats else None for stats in flood_buildings_stats]
#     buildings['flood_wd_75perc'] = [stats['75%'] if stats else None for stats in flood_buildings_stats]
#     buildings['flood_wd_max'] = [stats['max'] if stats else None for stats in flood_buildings_stats]

#     return buildings


# def compute_wd_summary(
#     buildings: gpd.GeoDataFrame,
#     provider: str,
#     include_stats: bool
# ) -> dict:
#     """
#     Compute summary statistics for flooded buildings.
#     This function will return a dictionary with aggregated statistics based on building type and class.
#     """
    
#     summary = dict()
    
#     def base_summary(gdf):
#         _base_summary = {
#             'total_buildings': len(gdf),
#             'flooded_buildings': int(gdf['is_flooded'].sum()),
#         }
#         if include_stats:
#             _base_summary.update({
#                 'flood_wd_min': float(np.nanmin(gdf['flood_wd_min'].values)),
#                 'flood_wd_25perc': float(np.nanpercentile(gdf['flood_wd_25perc'].values, 25)),
#                 'flood_wd_mean': float(np.nanmean(gdf['flood_wd_mean'].values)),
#                 'flood_wd_median': float(np.nanmedian(gdf['flood_wd_median'].values)),
#                 'flood_wd_75perc': float(np.nanpercentile(gdf['flood_wd_75perc'].values, 75)),
#                 'flood_wd_max': float(np.nanmax(gdf['flood_wd_max'].values))
#             })
#         _base_summary = {k: v if not np.isnan(v) else None for k, v in _base_summary.items()}
#         return _base_summary
    
#     summary['overall'] = base_summary(buildings)
    
#     class_column = None
#     if provider == 'OVERTURE':
#         class_column = 'subtype'
#     elif provider.startswith('RER-REST'):
#         class_column = 'rer_class'
#     else:
#         raise ValueError(f"Provider '{provider}' is not supported for summary computation. Available providers are: {_PROVIDERS}.")     # DOC: Should never happen, but just in case.
#     buildings[class_column] = buildings[class_column].fillna('other')
#     summary['classes'] = {
#         class_name: base_summary(class_gdf)
#         for class_name, class_gdf in buildings.groupby(class_column)
#     }
    
#     return summary


# DOC: Main function to compute flooded buildings

def compute_flood(
    waterdepth_filename: str,
    buildings_filename: str | None = None,
    wd_thresh: float = 0.5,
    bbox: tuple[float, float, float, float] | None = None,
    out: str | None = None,
    t_srs: str | None = None,
    provider: list[str] | None = None,
    feature_filters: dict[str, list[dict[str, list|str|int]]] | None = None,
    only_flood: bool = False,
    compute_stats: bool = False,
    compute_summary: bool = False
) -> str:
    
    """
    Main function to run the flooded buildings analysis.

    Parameters:
        waterdepth_filename (str): Path to the water depth raster file.
        buildings_filename (str | None): Path to the buildings vector file.
        bbox (tuple[float, float, float, float] | None): Bounding box for filtering.
        out (str | None): Output path for the results.
        t_srs (str | None): Target spatial reference system.
        providers (list[str] | None): List of data providers.
        feature_filters (dict[str, dict] | None): Filters for providers-features.

    Returns:
        str: Path to the output file containing flooded buildings.
    """

    # DOC: 0 — Print GECO logo
    module_logo.logo()
    
    
    # DOC: 1 — Validate args.
    Logger.info("# Validating input arguments ...")
    validated_args = module_args.validate_args(
        waterdepth_filename=waterdepth_filename,
        buildings_filename=buildings_filename,
        wd_thresh=wd_thresh,
        bbox=bbox,
        out=out,
        t_srs=t_srs,
        provider=provider,
        feature_filters=feature_filters,
        only_flood=only_flood,
        compute_stats=compute_stats,
        compute_summary=compute_summary
    )
    waterdepth_filename, buildings_filename, wd_thresh, bbox, out, t_srs, provider, feature_filters, only_flood, compute_stats, compute_summary = validated_args
    
    
    # DOC: 2 — Gather buildings
    Logger.info(f'# Gather buildings data ...')
    provider_buildings = module_retriever.retrieve_buildings(
        buildings_filename=buildings_filename,
        bbox=bbox,
        provider=provider
    )
    
    
    # DOC: 3 — Polygonize waterdepth (compute one time and reuse for multiple providers)
    Logger.info('# Processing water depth data (raster to polygon) ...')
    waterdepth_polygonized = _utils.polygonize_raster_valid_data(
        raster_filename=waterdepth_filename,
        band=1,
        mask_builder=lambda wd: wd > wd_thresh,    # Significant water depth threshold
        bbox=bbox
    )
    
    
    # DOC: 4 — Intersect buildings with water depth
    Logger.info('# Intersecting buildings with water depth ...')
    flooded_buildings = module_flood.get_flooded_buildings(
        waterdepth_mask = waterdepth_polygonized,
        buildings = provider_buildings
    )
        
    
    # DOC: 5 — Filter features
    Logger.info('# Filtering features ...')
    filtered_flooded_buildings = module_stats.filter_by_feature(
        gdf = flooded_buildings,
        feature_filters = feature_filters,
        only_flood = only_flood
    )
    Logger.info(f"## Filtered {len(filtered_flooded_buildings)} buildings out from {len(flooded_buildings)}.")
    
    
    # DOC: 6 — Compute water depth stats over flooded buildings
    if compute_stats:
        Logger.info('# Computing water depth stats over flooded buildings ...')
        filtered_flooded_buildings = module_stats.compute_wd_stats(
            waterdepth_filename=waterdepth_filename,
            waterdepth_mask=waterdepth_polygonized,
            buildings=filtered_flooded_buildings
        )
        Logger.info("## Water depth stats computed for flooded buildings.")
        
    
    # DOC: 7 — Compute summary if requested
    if compute_summary:
        Logger.info('# Computing summary statistics for flooded buildings ...')
        summary_stats = module_stats.compute_wd_summary(
            buildings=filtered_flooded_buildings,
            provider=provider,
            include_stats=compute_stats,
        )
        Logger.info("## Summary statistics computed for flooded buildings.") 
    
        
    
    # DOC: 8 — Return results
    Logger.info('# Preparing geojson output results ...')
    filtered_flooded_buildings = filtered_flooded_buildings.to_crs(t_srs)
    feature_collection = filtered_flooded_buildings.to_geo_dict()
    feature_collection['metadata'] = {
        'provider': provider,
        'buildings_count': len(filtered_flooded_buildings),
        'flooded_buildings_count': int(filtered_flooded_buildings['is_flooded'].sum()),
        'summary': summary_stats if compute_summary else None
    }
    feature_collection['crs'] = {
        "type": "name",
        "properties": {
            "name": f"urn:ogc:def:crs:{t_srs.replace(':', '::')}"  # REF: https://gist.github.com/sgillies/1233327 lines 256:271
        }
    }
    Logger.info(f"## Buildings feature collection prepared with {len(feature_collection['features'])} features.")
        
    
    # DOC: 9 — Save results to file
    Logger.info(f'# Saving results to {out} ...')
    if out.startswith('s3://'):
        out_tmp = f"{_utils.juststem(out)}.geojson"
        with open(out_tmp, 'w') as f:
            json.dump(feature_collection, f, indent=2)
        module_s3.s3_upload(filename = out_tmp,bucket = out)    
    else:
        with open(out, 'w') as f:
            json.dump(feature_collection, f, indent=2)
        
    
    Logger.info(f"## Results saved to {out}")
    
    return feature_collection



# DOC: Main function to run the flooded buildings analysis from command line

@click.command()
@click.option('--water', type=str, required=True, help='Path to the water depth raster file.')
@click.option('--building','--buildings', type=str, default=None, help='Path to the buildings vector file.')
@click.option('--wd_thresh', type=float, default=0.5, help='Water depth threshold for significant flooding (default: 0.5).')
@click.option('--bbox', type=float, nargs=4, default=None, help='Bounding box (minx, miny, maxx, maxy). If None, the total bounds of the water depth raster will be used.')
@click.option('--out', type=str, default=None, help='Output path for the results.')
@click.option('--t_srs', type=str, default=None, help='Target spatial reference system (EPSG code). If None, CRS of water depth raster will be used.')
@click.option('--provider', type=str, default=None, help='Building data provider (one of OVERTURE, RER-REST/*).')
@click.option('--filters', type=str, default=None, help='Filters for providers-features in JSON format.')
@click.option('--only_flood', is_flag=True, required=False, default=False, help="Only return flooded buildings (default: False).")
@click.option('--stats', is_flag=True, required=False, default=False, help="Compute water depth statistics for flooded buildings.")
@click.option('--summary', is_flag=True, required=False, default=False, help="Returns an additional metadata field with aggregated statistic based on building type and class. If true, stats will be computed as well.")

@click.option("--version", is_flag=True, required=False, default=False,
              help="Print version.")
@click.option("--debug", is_flag=True, required=False, default=False,
              help="Debug mode.")
@click.option("--verbose", required=False, is_flag=True, type=click.BOOL, default=False,
              help="verbose mode.")
def main(
    water,
    building,
    wd_thresh,
    bbox,
    out,
    t_srs,
    provider,
    filters,
    only_flood,
    stats,
    summary,
    version,
    debug,
    verbose
):
    """
    Main function to run the flooded buildings analysis from command line.
    
    Examples:
    1. safer-buildings --water tests\rimini-wd.tif --provider OVERTURE --filters "[{'subtype':'education', 'class': ['kindergarten','school']}, {'class':'parking'}]" --only_flood --stats
    2. safer-buildings --water tests\rimini-wd.tif --provider RER-REST/28/31/40 --filters "[{'ORDINE_NORMALIZZATO': ['Scuola primaria', 'Nido d\'infanzia']}, {'ISTITUZIONE_SCOLASTICA_RIF': 'IC ALIGHIERI'}]" --summary
    3. safer-buildings --water s3://saferplaces.co/Directed/process_out/SaferBuildingsService/rimini-wd.tif --buildings s3://saferplaces.co/Directed/process_out/SaferBuildingsService/Data/buildings-default-area__rer-rest_overture.geojson --out s3://saferplaces.co/Directed/process_out/SaferBuildingsService/rimini-wd-buildings.geojson --provider RER-REST/1 --summary
    
    In first example the water depth file is 'tests/rimini-wd.tif', the OVERTURE provider is used, and buildings are filtered is (subtype in ['education'] AND class in ['kindergarten', 'school']) OR (class in ['parking']).
    """
    
    Logger.info("# Starting flooded buildings analysis ...")
    time_start = time.time()
    
    if filters:
        filters = json.loads(filters.replace("'", '"'))  # Convert single quotes to double quotes for JSON parsing
        
    Logger.info("# CLI parameters:")
    Logger.info(f"## Water depth file: {water}")
    Logger.info(f"## Buildings file: {building}")
    Logger.info(f"## Water depth threshold: {wd_thresh}")
    Logger.info(f"## Bounding box: {bbox}")
    Logger.info(f"## Output file: {out}")
    Logger.info(f"## Target SRS: {t_srs}")
    Logger.info(f"## Provider: {provider}")
    Logger.info(f"## Feature filters: {filters}")
    Logger.info(f"## Only flood: {only_flood}")
    Logger.info(f"## Stats: {stats}")
    Logger.info(f"## Summary: {summary}")
    
    result = compute_flood(
        waterdepth_filename = water,
        buildings_filename = building,
        wd_thresh = wd_thresh,
        bbox = tuple(bbox) if bbox else None,
        out = out,
        t_srs = t_srs,
        provider = provider,
        feature_filters = filters,
        only_flood = only_flood,
        compute_stats = stats,
        compute_summary = summary
    )
    
    time_end = time.time()
    Logger.info(f"# Flooded buildings analysis completed in {time_end - time_start:.2f} seconds. Returned {len(result['features'])} features.")