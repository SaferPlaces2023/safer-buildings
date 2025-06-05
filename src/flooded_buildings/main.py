import os
import sys
import datetime
import json
import click
import time
from enum import Enum
import requests

import numpy as np
import pandas as pd

from osgeo import gdal, ogr, osr
from shapely.wkt import loads
from shapely.geometry import box, Point, Polygon
import geopandas as gpd

from eedem import downloadDEM as eedem_download

from collections.abc import MutableMapping

from . import utils



import logging
logging.basicConfig(format="%(asctime)s [%(levelname)s] %(message)s", stream=sys.stderr)
logger = logging.getLogger(__name__)    
logger.setLevel(logging.INFO)


# TODO: handle different crs (wd, bbox, overture, rest ecc...)
# TODO: buffer on proj crs (geo crs goes error !!)
# TODO: --list-providers arg shows detailed description of providers
# TODO: --wd-stats -> min,avg,max waterdepth around flooded buildings
# TODO: --summary -> add metadata with aggregate stats ( n buildings, n flooded buildings, min,avg,max wd for each category (provider dependend) of buildings)


# DOC: define costants


def _list_rest_layers():
    service_url = "https://servizigis.regione.emilia-romagna.it/geoags/rest/services/portale/saferplaces/MapServer"
    params = {
        "where": "1=1",
        "outFields": "*",
        "f": "json",
        "returnGeometry": "true"
    }
    headers = {
        "User-Agent": "QGIS",
        "Referer": "http://localhost"
    }
    response = requests.get(service_url, params=params, headers=headers)
    data = response.json()
    df_layers = pd.DataFrame(data['layers']).sort_values('id').reset_index(drop=True)
    df_layers['provider_name'] = df_layers['id'].apply(lambda service_id: f'REGIONE-EMILIA-ROMAGNA-{service_id}')
    return df_layers

RegioneEmiliaRomagnaLayers = _list_rest_layers()

_PROVIDERS = (
    'OVERTURE',
    * RegioneEmiliaRomagnaLayers.provider_name.to_list(),
)



# DOC: Aux functions


def validate_args(
    waterdepth_filename: str,
    buildings_filename: str | None = None,
    bbox: tuple[float, float, float, float] | None = None,
    out: str | None = None,
    t_srs: str | None = None,
    provider: list[str] | None = None,
    feature_filters: dict[str, dict] | None = None
):
    
    """
    Validate the input arguments for the flooded buildings analysis.
    """
        
    if waterdepth_filename is None:
        raise ValueError("waterdepth_filename must be provided.")
    if type(waterdepth_filename) is not str:
        raise TypeError("waterdepth_filename must be a string.")
    if os.path.isfile(waterdepth_filename) is False:
        raise FileNotFoundError(f"Water depth file not found: {waterdepth_filename}")
    
    if buildings_filename is not None:
        if type(buildings_filename) is not str:
            raise TypeError("buildings_filename must be a string.")
        if os.path.isfile(buildings_filename) is False:
            raise FileNotFoundError(f"Buildings file not found: {buildings_filename}")
        
    if bbox is not None:
        if len(bbox) != 4:
            raise ValueError("bbox must be a tuple of four floats (minx, miny, maxx, maxy).")
        if not all(isinstance(coord, (int, float)) for coord in bbox):
            raise TypeError("All coordinates in bbox must be int or float.")
    else:
        bbox = utils.get_raster_bounds(waterdepth_filename)
    
    if out is not None:
        if type(out) is not str:
            raise TypeError("out must be a string.") 
    else:
        out = os.path.join(os.getcwd(), f"flooded_buildings_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.geojson")
        
    if t_srs is not None:
        if type(t_srs) is not str:
            raise TypeError("t_srs must be a string.")
        if not t_srs.startswith("EPSG:"):
            raise ValueError("t_srs must be an EPSG code (e.g., 'EPSG:4326').")
        try:
            epsg_code = int(t_srs.split(":")[1])
        except ValueError:
            raise ValueError("Invalid EPSG code in t_srs.")
        srs = osr.SpatialReference()
        valid = srs.ImportFromEPSG(epsg_code)
        if valid != 0:
            raise ValueError(f"Invalid EPSG code: {epsg_code}")
    else:
        t_srs = utils.get_raster_crs(waterdepth_filename)
        
        
    if provider is None:
        if buildings_filename is not None:
            raise ValueError("A provider must be provided if buildings_filename is given.")
        else:
            provider = 'OVERTURE'
    if type(provider) is not str:
        raise TypeError("provider must be a string")
    if provider not in _PROVIDERS:
        raise ValueError(f"Invalid provider: {provider}. Valid providers are: {_PROVIDERS}.")
    
    
    if feature_filters is not None:
        if not isinstance(feature_filters, (dict, list)):
            raise TypeError("feature_filters must be a dictionary or a list.")
        if isinstance(feature_filters, dict):
            feature_filters = [feature_filters]
        for idxf, filters in enumerate(feature_filters):
            for f_key, f_value in filters.items():
                if not isinstance(f_key, str):
                    raise TypeError("Filter keys must be strings.")
                if not isinstance(f_value, (list, str, int)):
                    raise TypeError("Filter values must be a list, string, or integer.")
                if not isinstance(f_value, list):
                    f_value = [f_value]
                filters[f_key] = f_value
            feature_filters[idxf] = filters                            
    else:
        feature_filters = []
        
    print("## Input arguments validated successfully.")
    print(f"### Water depth file: {waterdepth_filename}")
    print(f"### Buildings file: {buildings_filename}")
    print(f"### Bounding box: {bbox}")
    print(f"### Output file: {out}")
    print(f"### Target SRS: {t_srs}")
    print(f"### Providers: {provider}")
    print(f"### Feature filters: {feature_filters}")
        
    return waterdepth_filename, buildings_filename, bbox, out, t_srs, provider, feature_filters


def retrieve_buildings(
    buildings_filename: str | None, 
    bbox: tuple[float, float, float, float] | None,
    provider: str | None
) -> gpd.GeoDataFrame:
    
    """
    Retrieve buildings data from specified providers.
    If buildings_filename is provided, it will be used directly.
    If not, it will download buildings data from the specified providers.
    """
    
    if buildings_filename is not None:
        provider_buildings = gpd.read_file(buildings_filename)
        print(f"## Using provided buildings data from {buildings_filename}. Found {len(provider_buildings)} buildings.")
    
    else:
        print(f"## Retrieving buildings data from provider: {provider} ...")
            
        buildings_filename = utils.temp_filename(ext='shp', prefix=f'{provider}_buildings')
        
        if provider == 'OVERTURE':
            _ = eedem_download(
                dataset = 'OVERTURE/BUILDINGS',
                bbox = utils.shapely_bbox_2_eedem_bbox(box(*bbox)),
                band=None,
                out = buildings_filename,
                dmg = True
            )
            provider_buildings = gpd.read_file(buildings_filename)
            
        elif provider.startswith('REGIONE-EMILIA-ROMAGNA'):
            # DOC: Retrieve from all subservices of the requested one.
            service_id = int(provider.split('-')[-1])
            service_ids = [service_id]
            while any([RegioneEmiliaRomagnaLayers[RegioneEmiliaRomagnaLayers.id == service_id].iloc[0].subLayerIds is not None for service_id in service_ids]):
                for service_id in service_ids:
                    sub_layers = RegioneEmiliaRomagnaLayers[RegioneEmiliaRomagnaLayers.id == service_id].iloc[0].subLayerIds
                    if sub_layers is not None:
                        service_ids.remove(service_id)
                        service_ids.extend(sub_layers)
                        
            def rest_service_retrieve(service_id):
                url = f"https://servizigis.regione.emilia-romagna.it/geoags/rest/services/portale/saferplaces/MapServer/{service_id}/query"
                params = {
                    "geometry": ','.join([str(b) for b in bbox]),
                    "geometryType": "esriGeometryEnvelope",
                    "inSR": "4326",
                    "outSR": "7791",
                    "spatialRel": "esriSpatialRelIntersects",
                    "where": "1=1",
                    "outFields": "*",
                    "f": "json",
                    "returnGeometry": "true"
                }
                headers = {
                    "User-Agent": "QGIS",
                    "Referer": "http://localhost"
                }
                response = requests.get(url, params=params, headers=headers)
                data = response.json()
                
                features =[f['attributes'] for f in data['features']]
                geometries = [Point(f['geometry']['x'], f['geometry']['y']) for f in data['features']]
                
                buffer_meters = 20

                rest_gdf = gpd.GeoDataFrame(features, geometry=geometries, crs=f"EPSG:{data['spatialReference']['wkid']}")
                rest_gdf['geometry'] = rest_gdf.buffer(buffer_meters)
                rest_gdf = rest_gdf.to_crs(epsg=4326)
                return rest_gdf
            
            provider_buildings = pd.concat([rest_service_retrieve(service_id) for service_id in service_ids], ignore_index=True)
            provider_buildings.to_file(buildings_filename, driver='ESRI Shapefile', index=False)
        
        print(f"### Buildings data retrieved from {provider} saved at {buildings_filename}. (Found {len(provider_buildings)} buildings)")        
        
    return provider_buildings


def get_waterdepth_mask(
    waterdepth_filename: str,
    mask_builder: callable = lambda wd: wd > 0.5,  # Default mask for significant water depth
) -> gpd.GeoDataFrame:
    
    """
    Create a mask from the water depth raster file.
    The mask is created by polygonizing the raster data based on the provided mask_builder function.
    """
    
    waterdepth_polygonized = utils.polygonize_raster_valid_data(
        raster_filename=waterdepth_filename,
        band=1,
        mask_builder=mask_builder
    )
    
    return waterdepth_polygonized


def get_flooded_buildings(
    waterdepth_filename: str,
    waterdepth_mask: gpd.GeoDataFrame | None,
    buildings: gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:
    
    """ 
    Get flooded buildings by intersecting water depth polygons with buildings.
    """
    
    buildings = utils.ensure_geodataframe_crs(buildings, utils.get_raster_crs(waterdepth_filename))
    
    # TODO: crop to bbox if defined (assume user privided bbox is 4326)
    
    if waterdepth_mask is None:
        waterdepth_mask = get_waterdepth_mask(
            waterdepth_filename=waterdepth_filename,
            mask_builder=lambda wd: wd > 0.5  # Significant water depth threshold
        )
    
    # Intersect with water depth polygons
    buildings['__tmp_identifier__'] = buildings.index.to_list()
    flooded_buildings = gpd.overlay(buildings, waterdepth_mask, how='intersection')
    buildings['is_flooded'] = buildings['__tmp_identifier__'].apply(lambda tmp_id: tmp_id in flooded_buildings['__tmp_identifier__'].to_list())
    buildings.drop(columns=['__tmp_identifier__'], inplace=True)
    
    print(f"### Found {len(flooded_buildings)} flooded buildings out of {len(buildings)} total buildings.")
    
    return buildings


def filter_by_feature(
    gdf: gpd.GeoDataFrame,
    feature_filters: list[dict[str, list]]
):
    """
    Filter a GeoDataFrame by a list of filters.
    """
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
    


# DOC: Main function to compute flooded buildings

def compute_flood(
    waterdepth_filename: str,
    buildings_filename: str | None = None,
    bbox: tuple[float, float, float, float] | None = None,
    out: str | None = None,
    t_srs: str | None = None,
    provider: list[str] | None = None,
    feature_filters: dict[str, list[dict[str, list|str|int]]] | None = None,
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
    
    
    # DOC: 1 — Validate args.
    print("# Validating input arguments ...")
    waterdepth_filename, buildings_filename, bbox, out, t_srs, provider, feature_filters = validate_args(
        waterdepth_filename=waterdepth_filename,
        buildings_filename=buildings_filename,
        bbox=bbox,
        out=out,
        t_srs=t_srs,
        provider=provider,
        feature_filters=feature_filters
    )
    
    
    # DOC: 2 — Retrieve buildings
    print(f'# Retrieving buildings data from provider {provider} ...')
    provider_buildings = retrieve_buildings(
        buildings_filename=buildings_filename,
        bbox=bbox,
        provider=provider
    )
    
    
    # DOC: 3 — Polygonize waterdepth (compute one time and reuse for multiple providers)
    print('# Processing water depth data (raster to polygon) ...')
    waterdepth_polygonized = utils.polygonize_raster_valid_data(
        raster_filename=waterdepth_filename,
        band=1,
        mask_builder=lambda wd: wd > 0.5    # Significant water depth threshold
    )
    
    # TODO: crop to bbox if defined (assume user provided bbox is 4326)
    
    
    # DOC: 4 — Intersect buildings with water depth
    print('# Intersecting buildings with water depth ...')
    provider_flooding = get_flooded_buildings(
        waterdepth_filename = waterdepth_filename,
        waterdepth_mask = waterdepth_polygonized,
        buildings = provider_buildings
    )
        
    
    # DOC: 5 — Filter features
    print('# Filtering features ...')
    filtered_provider_gdf = filter_by_feature(
        gdf = provider_flooding,
        feature_filters = feature_filters
    )
    print(f"## Filtered {len(filtered_provider_gdf)} buildings out from {len(provider_flooding)}.")
        
    
    # DOC: 6 — Return results
    print('# Preparing geojson output results ...')
    provider_feature_collection = filtered_provider_gdf.to_geo_dict()
    provider_feature_collection['metadata'] = {
        'provider': provider,
        'buildings_count': len(filtered_provider_gdf),
        'flooded_buildings_count': len(filtered_provider_gdf['is_flooded'])
    }
    provider_feature_collection['crs'] = {
        "type": "name",
        "properties": {
            "name": f"urn:ogc:def:crs:{t_srs.replace(':', '::')}"  # REF: https://gist.github.com/sgillies/1233327 lines 256:271
        }
    }
    print(f"## Buildings feature collection prepared with {len(provider_feature_collection['features'])} features.")
        
    
    # DOC: 7 — Save results to file
    print(f'# Saving results to {out} ...')
    out_provider_fname = f'{out}__{provider}.geojson'
    with open(out_provider_fname, 'w') as f:
        json.dump(provider_feature_collection, f, indent=2)
    print(f"## Results saved to {out_provider_fname}")
    
    return provider_feature_collection



# DOC: Main function to run the flooded buildings analysis from command line

@click.command()
@click.option('--wd', type=click.Path(exists=True), required=True, help='Path to the water depth raster file.')
@click.option('--buildings', type=click.Path(exists=True), default=None, help='Path to the buildings vector file.')
@click.option('--bbox', type=float, nargs=4, default=None, help='Bounding box (minx, miny, maxx, maxy).')
@click.option('--out', type=click.Path(), default=None, help='Output path for the results.')
@click.option('--t_srs', type=str, default=None, help='Target spatial reference system (EPSG code).')
@click.option('--provider', type=str, default=None, help='Building data provider (one of OVERTURE, REGIONE-EMILIA-ROMAGNA-*).')
@click.option('--filters', type=str, default=None, help='Filters for providers-features in JSON format.')
def main(wd, buildings, bbox, out, t_srs, provider, filters):
    """
    Main function to run the flooded buildings analysis from command line.
    
    Examples:
    1. safer-buildings --wd tests\rimini-wd.tif --provider OVERTURE --filters "[{'subtype':'education', 'class': ['kindergarten','school']}, {'class':'parking'}]"
    2. safer-buildings --wd tests\rimini-wd.tif --provider REGIONE-EMILIA-ROMAGNA-30 --filters "[{'ORDINE_NORMALIZZATO': ['Scuola primaria', 'Nido d\'infanzia']}, {'ISTITUZIONE_SCOLASTICA_RIF': 'IC ALIGHIERI'}]"
    
    In first example the water depth file is 'tests/rimini-wd.tif', the OVERTURE provider is used, and buildings are filtered is (subtype in ['education'] AND class in ['kindergarten', 'school']) OR (class in ['parking']).
    """
    
    print("# Starting flooded buildings analysis ...")
    time_start = time.time()
    
    if filters:
        filters = json.loads(filters.replace("'", '"'))  # Convert single quotes to double quotes for JSON parsing
        
    print("# CLI parameters:")
    print(f"## Water depth file: {wd}")
    print(f"## Buildings file: {buildings}")
    print(f"## Bounding box: {bbox}")
    print(f"## Output file: {out}")
    print(f"## Target SRS: {t_srs}")
    print(f"## Provider: {provider}")
    print(f"## Feature filters: {filters}")
    
    result = compute_flood(
        waterdepth_filename = wd,
        buildings_filename = buildings,
        bbox = tuple(bbox) if bbox else None,
        out = out,
        t_srs = t_srs,
        provider = provider,
        feature_filters = filters
    )
    
    time_end = time.time()
    print(f"# Flooded buildings analysis completed in {time_end - time_start:.2f} seconds.")