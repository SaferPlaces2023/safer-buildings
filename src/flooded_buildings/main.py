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

from osgeo import gdal, ogr, osr
from shapely.wkt import loads
from shapely.geometry import box, Point, Polygon, MultiPolygon, LineString, MultiLineString
import geopandas as gpd

from eedem import downloadDEM as eedem_download

from . import utils



logging.basicConfig(format="%(asctime)s [%(levelname)s] %(message)s", stream=sys.stderr)
logger = logging.getLogger(__name__)    
logger.setLevel(logging.INFO)



# TODO: --list-providers → shows detailed description of providers
# TODO: --summary → add metadata with aggregate stats ( n buildings, n flooded buildings, min,avg,max wd for each category (provider dependend) of buildings)



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
    wd_thresh: float = 0.5,
    bbox: tuple[float, float, float, float] | None = None,
    out: str | None = None,
    t_srs: str | None = None,
    provider: list[str] | None = None,
    feature_filters: dict[str, dict] | None = None,
    only_flood: bool = False,
    compute_stats: bool = False
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
        
    if wd_thresh is None:
        wd_thresh = 0.5
    if type(wd_thresh) not in (int, float):
        raise TypeError("wd_thresh must be a float or int.")
    if wd_thresh < 0:
        raise ValueError("wd_thresh must be a non-negative float or int.")
        
    if bbox is not None:
        if len(bbox) != 4:
            raise ValueError("bbox must be a tuple of four floats (minx, miny, maxx, maxy).")
        if not all(isinstance(coord, (int, float)) for coord in bbox):
            raise TypeError("All coordinates in bbox must be int or float.")
        bbox = gpd.GeoDataFrame({'geometry': [box(*bbox)]}, crs="EPSG:4326")
    else:
        bbox = gpd.GeoDataFrame({'geometry': [box(*utils.get_raster_bounds(waterdepth_filename))]}, crs=utils.get_raster_crs(waterdepth_filename))
    
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
        
    if only_flood is None:
        only_flood = False
    if type(only_flood) is not bool:
        raise TypeError("only_flood must be a boolean value.")
        
    if compute_stats is None:
        compute_stats = False
    if type(compute_stats) is not bool:
        raise TypeError("compute_stats must be a boolean value.")
        
    print("## Input arguments validated successfully.")
    print(f"### Water depth file: {waterdepth_filename}")
    print(f"### Buildings file: {buildings_filename}")
    print(f"### Water depth threshold: {wd_thresh}")
    print(f"### Bounding box (total bounds): {bbox.total_bounds}")
    print(f"### Output file: {out}")
    print(f"### Target SRS: {t_srs}")
    print(f"### Providers: {provider}")
    print(f"### Feature filters: {feature_filters}")
    print(f"### Only flood: {only_flood}")
    print(f"### Compute stats: {compute_stats}")
        
    return waterdepth_filename, buildings_filename, wd_thresh, bbox, out, t_srs, provider, feature_filters, only_flood, compute_stats


def retrieve_buildings(
    buildings_filename: str | None, 
    bbox: tuple[float, float, float, float],
    provider: str
) -> gpd.GeoDataFrame:
    
    """
    Retrieve buildings data from specified providers.
    If buildings_filename is provided, it will be used directly.
    If not, it will download buildings data from the specified providers.
    """
    
    if buildings_filename is not None:
        provider_buildings = gpd.read_file(buildings_filename)
        bbox = gpd.GeoDataFrame({'geometry': [box(*bbox)]}, crs="EPSG:4326").to_crs(utils.get_geodataframe_crs(provider_buildings))
        provider_buildings = gpd.overlay(provider_buildings, bbox, how='intersection')
        print(f"## Using provided buildings data from {buildings_filename}. Found {len(provider_buildings)} buildings.")
    
    else:
        print(f"## Retrieving buildings data from provider: {provider} ...")
            
        buildings_filename = utils.temp_filename(ext='shp', prefix=f'{provider}_buildings')
        
        if provider == 'OVERTURE':
            bbox4326 = bbox.to_crs(epsg=4326).geometry.iloc[0]
            print(utils.shapely_bbox_2_eedem_bbox(bbox4326),)
            _ = eedem_download(
                dataset = 'OVERTURE/BUILDINGS',
                bbox = utils.shapely_bbox_2_eedem_bbox(bbox4326),
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
                bbox4326 = bbox.to_crs("EPSG:4326").geometry.total_bounds
                params = {
                    "geometry": ','.join([str(b) for b in bbox4326]),
                    "geometryType": "esriGeometryEnvelope",
                    "inSR": "4326",
                    "outSR": "4326",
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
                
                do_buffer = False
                geometry_type = RegioneEmiliaRomagnaLayers[RegioneEmiliaRomagnaLayers.id == service_id].iloc[0].geometryType
                if geometry_type == 'esriGeometryPoint':
                    geometries = [
                        Point(f['geometry']['x'], f['geometry']['y']) 
                        for f in data['features']
                    ]
                    do_buffer = True
                elif geometry_type == 'esriGeometryPolygon':
                    geometries = [
                        MultiPolygon(
                            [Polygon([point for point in ring ])
                            for ring in f['geometry']['rings']]
                        ) 
                        for f in data['features']
                    ]
                elif geometry_type == 'esriGeometryPolyline':
                    geometries = [
                        MultiLineString(
                            [LineString([point for point in path ])
                            for path in f['geometry']['paths']]
                        ) 
                        for f in data['features']
                    ]
                    do_buffer = True
                else:
                    raise ValueError(f"Unsupported geometry type for service {service_id}: {RegioneEmiliaRomagnaLayers[RegioneEmiliaRomagnaLayers.id == service_id].iloc[0].geometryType}")

                rest_gdf = gpd.GeoDataFrame(features, geometry=geometries, crs=f"EPSG:{data['spatialReference']['wkid']}")
                rest_gdf['service_id'] = service_id
                rest_gdf['service_class'] = RegioneEmiliaRomagnaLayers[RegioneEmiliaRomagnaLayers.id == service_id].name.iloc[0]
                
                if do_buffer:
                    buffer_meters = 20
                    rest_gdf = rest_gdf.to_crs("EPSG:7791")
                    rest_gdf['geometry'] = rest_gdf.buffer(buffer_meters)
                    
                return rest_gdf
        
            provider_buildings = pd.concat([rest_service_retrieve(service_id) for service_id in service_ids], ignore_index=True)
            provider_buildings.to_file(buildings_filename, driver='ESRI Shapefile', index=False)
        
        else:
            raise ValueError(f"Provider '{provider}' is not supported. Available providers are: {_PROVIDERS}.")

        
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
    waterdepth_mask: gpd.GeoDataFrame | None,
    buildings: gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:
    
    """ 
    Get flooded buildings by intersecting water depth polygons with buildings.
    """
    
    buildings = utils.ensure_geodataframe_crs(buildings, utils.get_geodataframe_crs(waterdepth_mask))
    
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
    

def compute_wd_stats(
    waterdepth_filename: str,
    waterdepth_mask: gpd.GeoDataFrame,
    buildings: gpd.GeoDataFrame
):
    """
    Compute statistics on water depth around flooded buildings.
    """
    
    def building_wd_stats(building):
        flood_area = gpd.overlay(waterdepth_mask, gpd.GeoDataFrame({'geometry': [building.geometry]}), how='intersection') if building.is_flooded else None
        if flood_area is None or flood_area.empty:
            return None
        else:
            flood_area_values = utils.raster_sample_area(waterdepth_filename, flood_area.geometry.iloc[0])
            flood_area_stats = dict(pd.Series(flood_area_values).describe())
            return flood_area_stats
               
    flood_buildings_stats = [building_wd_stats(building) for _, building in buildings.iterrows()]
    buildings['flood_wd_min'] = [stats['min'] if stats else None for stats in flood_buildings_stats]
    buildings['flood_wd_25perc'] = [stats['25%'] if stats else None for stats in flood_buildings_stats]
    buildings['flood_wd_mean'] = [stats['mean'] if stats else None for stats in flood_buildings_stats]
    buildings['flood_wd_75perc'] = [stats['75%'] if stats else None for stats in flood_buildings_stats]
    buildings['flood_wd_max'] = [stats['max'] if stats else None for stats in flood_buildings_stats]
    
    return buildings



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
    compute_stats: bool = False
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
    waterdepth_filename, buildings_filename, wd_thresh, bbox, out, t_srs, provider, feature_filters, only_flood, compute_stats = validate_args(
        waterdepth_filename=waterdepth_filename,
        buildings_filename=buildings_filename,
        wd_thresh=wd_thresh,
        bbox=bbox,
        out=out,
        t_srs=t_srs,
        provider=provider,
        feature_filters=feature_filters,
        only_flood=only_flood,
        compute_stats=compute_stats
    )
    
    
    # DOC: 2 — Gather buildings
    print(f'# Gather buildings data ...')
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
        mask_builder=lambda wd: wd > wd_thresh,    # Significant water depth threshold
        bbox=bbox
    )
    
    
    # DOC: 4 — Intersect buildings with water depth
    print('# Intersecting buildings with water depth ...')
    flooded_buildings = get_flooded_buildings(
        waterdepth_mask = waterdepth_polygonized,
        buildings = provider_buildings
    )
        
    
    # DOC: 5 — Filter features
    print('# Filtering features ...')
    filtered_flooded_buildings = filter_by_feature(
        gdf = flooded_buildings,
        feature_filters = feature_filters
    )
    print(f"## Filtered {len(filtered_flooded_buildings)} buildings out from {len(flooded_buildings)}.")
    
    # DOC:: 5.1 — Filter only flooded buildings if requested
    if only_flood:
        print('# Filtering only flooded buildings ...')
        filtered_flooded_buildings = filtered_flooded_buildings[filtered_flooded_buildings['is_flooded']]
        print(f"## Only flooded buildings retained: {len(filtered_flooded_buildings)} out of {len(flooded_buildings)}.")
    
    
    # DOC: 6 — Compute water depth stats over flooded buildings
    if compute_stats:
        print('# Computing water depth stats over flooded buildings ...')
        filtered_flooded_buildings = compute_wd_stats(
            waterdepth_filename=waterdepth_filename,
            waterdepth_mask=waterdepth_polygonized,
            buildings=filtered_flooded_buildings
        )
        print("## Water depth stats computed for flooded buildings.")
        
    
    # DOC: 7 — Return results
    print('# Preparing geojson output results ...')
    filtered_flooded_buildings = filtered_flooded_buildings.to_crs(t_srs)
    feature_collection = filtered_flooded_buildings.to_geo_dict()
    feature_collection['metadata'] = {
        'provider': provider,
        'buildings_count': len(filtered_flooded_buildings),
        'flooded_buildings_count': len(filtered_flooded_buildings['is_flooded'])
    }
    feature_collection['crs'] = {
        "type": "name",
        "properties": {
            "name": f"urn:ogc:def:crs:{t_srs.replace(':', '::')}"  # REF: https://gist.github.com/sgillies/1233327 lines 256:271
        }
    }
    print(f"## Buildings feature collection prepared with {len(feature_collection['features'])} features.")
        
    
    # DOC: 8 — Save results to file
    print(f'# Saving results to {out} ...')
    out_provider_fname = f'{out}__{provider}.geojson'
    with open(out_provider_fname, 'w') as f:
        json.dump(feature_collection, f, indent=2)
    print(f"## Results saved to {out_provider_fname}")
    
    return feature_collection



# DOC: Main function to run the flooded buildings analysis from command line

@click.command()
@click.option('--wd', type=click.Path(exists=True), required=True, help='Path to the water depth raster file.')
@click.option('--buildings', type=click.Path(exists=True), default=None, help='Path to the buildings vector file.')
@click.option('--wd-thresh', type=float, default=0.5, help='Water depth threshold for significant flooding (default: 0.5).')
@click.option('--bbox', type=float, nargs=4, default=None, help='Bounding box (minx, miny, maxx, maxy).')
@click.option('--out', type=click.Path(), default=None, help='Output path for the results.')
@click.option('--t_srs', type=str, default=None, help='Target spatial reference system (EPSG code).')
@click.option('--provider', type=str, default=None, help='Building data provider (one of OVERTURE, REGIONE-EMILIA-ROMAGNA-*).')
@click.option('--filters', type=str, default=None, help='Filters for providers-features in JSON format.')
@click.option('--only_flood', is_flag=True, required=False, default=False, help="Only return flooded buildings (default: False).")
@click.option('--stats', is_flag=True, required=False, default=False, help="Compute water depth statistics for flooded buildings.")
def main(wd, buildings, wd_thresh, bbox, out, t_srs, provider, filters, only_flood, stats):
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
    print(f"## Water depth threshold: {wd_thresh}")
    print(f"## Bounding box: {bbox}")
    print(f"## Output file: {out}")
    print(f"## Target SRS: {t_srs}")
    print(f"## Provider: {provider}")
    print(f"## Feature filters: {filters}")
    print(f"## Only flood: {only_flood}")
    print(f"## Stats: {stats}")
    
    result = compute_flood(
        waterdepth_filename = wd,
        buildings_filename = buildings,
        wd_thresh = wd_thresh,
        bbox = tuple(bbox) if bbox else None,
        out = out,
        t_srs = t_srs,
        provider = provider,
        feature_filters = filters,
        only_flood = only_flood,
        compute_stats = stats
    )
    
    time_end = time.time()
    print(f"# Flooded buildings analysis completed in {time_end - time_start:.2f} seconds. Returned {len(result['features'])} features.")