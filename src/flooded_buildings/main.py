import os
import sys
import datetime
import json
import click
import time
from enum import Enum

import numpy as np
import pandas as pd

from osgeo import gdal, ogr, osr
from shapely.wkt import loads
from shapely.geometry import box
import geopandas as gpd

from eedem import downloadDEM as eedem_download

from . import utils


import logging
logging.basicConfig(format="%(asctime)s [%(levelname)s] %(message)s", stream=sys.stderr)
logger = logging.getLogger(__name__)    
logger.setLevel(logging.INFO)



# DOC: define costants


class Providers(str, Enum):
    OVERTURE = 'OVERTURE'
    REGIONE_EMILIA_ROMAGNA = 'REGIONE_EMILIA_ROMAGNA'
    
class OvertureSubtypes(str, Enum):
    agricultural = 'agricultural'
    civic = 'civic'
    commercial = 'commercial'
    education = 'education'
    entertainment = 'entertainment'
    industrial = 'industrial'
    medical = 'medical'
    military = 'military'
    outbuilding = 'outbuilding'
    religious = 'religious'
    residential = 'residential'
    service = 'service'
    transportation = 'transportation'
    
class OvertureClass(str, Enum):
    agricultural = 'agricultural'
    allotment_house = 'allotment_house'
    apartments = 'apartments'
    barn = 'barn'
    beach_hut = 'beach_hut'
    boathouse = 'boathouse'
    bridge_structure = 'bridge_structure'
    bungalow = 'bungalow'
    bunker = 'bunker'
    cabin = 'cabin'
    carport = 'carport'
    cathedral = 'cathedral'
    chapel = 'chapel'
    church = 'church'
    civic = 'civic'
    college = 'college'
    commercial = 'commercial'
    cowshed = 'cowshed'
    detached = 'detached'
    digester = 'digester'
    dormitory = 'dormitory'
    dwelling_house = 'dwelling_house'
    factory = 'factory'
    farm = 'farm'
    farm_auxiliary = 'farm_auxiliary'
    fire_station = 'fire_station'
    garage = 'garage'
    garages = 'garages'
    ger = 'ger'
    glasshouse = 'glasshouse'
    government = 'government'
    grandstand = 'grandstand'
    greenhouse = 'greenhouse'
    guardhouse = 'guardhouse'
    hangar = 'hangar'
    hospital = 'hospital'
    hotel = 'hotel'
    house = 'house'
    houseboat = 'houseboat'
    hut = 'hut'
    industrial = 'industrial'
    kindergarten = 'kindergarten'
    kiosk = 'kiosk'
    library = 'library'
    manufacture = 'manufacture'
    military = 'military'
    monastery = 'monastery'
    mosque = 'mosque'
    office = 'office'
    outbuilding = 'outbuilding'
    parking = 'parking'
    pavilion = 'pavilion'
    post_office = 'post_office'
    presbytery = 'presbytery'
    public = 'public'
    religious = 'religious'
    residential = 'residential'
    retail = 'retail'
    roof = 'roof'
    school = 'school'
    semi = 'semi'
    semidetached_house = 'semidetached_house'
    service = 'service'
    shed = 'shed'
    shrine = 'shrine'
    silo = 'silo'
    slurry_tank = 'slurry_tank'
    sports_centre = 'sports_centre'
    sports_hall = 'sports_hall'
    stable = 'stable'
    stadium = 'stadium'
    static_caravan = 'static_caravan'
    stilt_house = 'stilt_house'
    storage_tank = 'storage_tank'
    sty = 'sty'
    supermarket = 'supermarket'
    synagogue = 'synagogue'
    temple = 'temple'
    terrace = 'terrace'
    toilets = 'toilets'
    train_station = 'train_station'
    transformer_tower = 'transformer_tower'
    transportation = 'transportation'
    trullo = 'trullo'
    university = 'university'
    warehouse = 'warehouse'
    wayside_shrine = 'wayside_shrine'



# DOC: Aux functions


def validate_args(
    waterdepth_filename: str,
    buildings_filename: str | None = None,
    bbox: tuple[float, float, float, float] | None = None,
    out: str | None = None,
    t_srs: str | None = None,
    providers: list[str] | None = None,
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
        
    if providers is not None:
        if not isinstance(providers, (list, tuple)):
            raise TypeError("providers must be a list of strings.")
        providers = list(providers)  # Ensure providers is a list
        if not all(isinstance(provider, str) for provider in providers):
            raise TypeError("All providers must be strings.")
        if any(provider not in Providers._member_names_ for provider in providers):
            raise ValueError(f"Invalid provider(s) in providers: {', '.join(providers)}. Valid providers are: {Providers._member_names_}.")
    else:
        providers = Providers._member_names_
    
    if feature_filters is not None:
        if not isinstance(feature_filters, dict):
            raise TypeError("feature_filters must be a dictionary.")
        for provider, provider_filters in feature_filters.items():
            if not isinstance(provider, str):
                raise TypeError("Provider names in feature_filters must be strings.")
            if provider not in providers:
                raise ValueError(f"Provider '{provider}' in feature_filters is not in the list of providers: {providers}.")
            if provider not in Providers._member_names_:
                raise ValueError(f"Invalid provider '{provider}' in feature_filters. Valid providers are: {Providers._member_names_}.")
            if not isinstance(provider_filters, (dict, list)):
                raise TypeError("Filters for each provider must be a dictionary or a list.")
            if isinstance(provider_filters, dict):
                provider_filters = [provider_filters]
            
            for idxf, filters in enumerate(provider_filters):
                for f_key, f_value in filters.items():
                    if not isinstance(f_key, str):
                        raise TypeError("Filter keys must be strings.")
                    if not isinstance(f_value, (list, str, int)):
                        raise TypeError("Filter values must be a list, string, or integer.")
                    if not isinstance(f_value, list):
                        f_value = [f_value]
                    if provider == Providers.OVERTURE:
                        if f_key == 'subtype':
                            if any(subtype not in OvertureSubtypes._member_names_ for subtype in f_value):
                                raise ValueError(f"Invalid subtype(s) in feature_filters for OVERTURE: {', '.join(f_value)}. Valid subtypes are: {OvertureSubtypes._member_names_}.")
                        elif f_key == 'class':
                            if any(cl not in OvertureClass._member_names_ for cl in f_value):
                                raise ValueError(f"Invalid class(es) in feature_filters for OVERTURE: {', '.join(f_value)}. Valid classes are: {OvertureClass._member_names_}.")
                        else:
                            raise ValueError(f"Invalid filter key '{f_key}' for provider '{provider}'. Valid keys are: 'subtype', 'class'.")
                    # TODO: Implement other providers
                    
                    filters[f_key] = f_value
                provider_filters[idxf] = filters
            feature_filters[provider] = provider_filters
                            
    else:
        feature_filters = {provider: [] for provider in providers}
        
    print("## Input arguments validated successfully.")
    print(f"### Water depth file: {waterdepth_filename}")
    print(f"### Buildings file: {buildings_filename}")
    print(f"### Bounding box: {bbox}")
    print(f"### Output file: {out}")
    print(f"### Target SRS: {t_srs}")
    print(f"### Providers: {providers}")
    print(f"### Feature filters: {feature_filters}")
        
    return waterdepth_filename, buildings_filename, bbox, out, t_srs, providers, feature_filters


def retrieve_buildings(
    buildings_filename: str | None, 
    bbox: tuple[float, float, float, float] | None,
    providers: dict[str, str] | None
) -> str:
    
    """
    Retrieve buildings data from specified providers.
    If buildings_filename is provided, it will be used directly.
    If not, it will download buildings data from the specified providers.
    """
    
    if buildings_filename is not None:
        buildings = gpd.read_file(buildings_filename)
        return buildings
    
    else:
        provider_buildings = dict() 
        
        for provider in providers:
            
            print(f"## Retrieving buildings data from provider: {provider} ...")
            
            if provider == Providers.OVERTURE:
                buildings_filename = utils.temp_filename(ext='shp', prefix='overture_buildings_')
                _ = eedem_download(
                    dataset = 'OVERTURE/BUILDINGS',
                    bbox = utils.shapely_bbox_2_eedem_bbox(box(*bbox)),
                    band=None,
                    out = buildings_filename,
                    dmg = True
                )
                provider_buildings[Providers.OVERTURE.name] = buildings_filename
                
            # TODO: Implement other providers
            
            print(f"### Buildings data retrieved from {provider} saved at {buildings_filename}. (Found {len(gpd.read_file(buildings_filename))} buildings)")
            provider_buildings[provider] = buildings_filename
        
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
    buildings_filename: str | gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:
    
    """ 
    Get flooded buildings by intersecting water depth polygons with buildings.
    """
    
    buildings = gpd.read_file(buildings_filename) if type(buildings_filename) is str else buildings_filename
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
    
    print(f"### Found {len(flooded_buildings)} flooded buildings out of {len(buildings)} total buildings in {buildings_filename}.")
    
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
    providers: list[str] | None = None,
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
    waterdepth_filename, buildings_filename, bbox, out, t_srs, providers, feature_filters = validate_args(
        waterdepth_filename=waterdepth_filename,
        buildings_filename=buildings_filename,
        bbox=bbox,
        out=out,
        t_srs=t_srs,
        providers=providers,
        feature_filters=feature_filters
    )
    
    
    # DOC: 2 — Retrieve buildings
    print('# Retrieving buildings data ...')
    buildings_providers_filenames = retrieve_buildings(
        buildings_filename=buildings_filename,
        bbox=bbox,
        providers=providers
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
    floodings_providers_gdfs = dict()
    for provider, buildings_filename in buildings_providers_filenames.items():
        print(f"## Processing intersection with provider: {provider} ...")
        floodings_providers_gdfs[provider] = get_flooded_buildings(
            waterdepth_filename = waterdepth_filename,
            waterdepth_mask = waterdepth_polygonized,
            buildings_filename = buildings_filename
        )
        
    
    # DOC: 5 — Filter features
    print('# Filtering features ...')
    for provider, provider_filters in feature_filters.items():
        print(f"## Filtering provider: {provider} with filters: {provider_filters}")
        provider_gdf = floodings_providers_gdfs[provider]
        filtered_provider_gdf = filter_by_feature(
            gdf = provider_gdf,
            feature_filters = provider_filters
        )
        print(f"### Filtered {len(filtered_provider_gdf)} buildings out from {len(provider_gdf)} for provider {provider}.")
        floodings_providers_gdfs[provider] = filtered_provider_gdf
        
    
    # DOC: 6 — Return results
    print('# Preparing geojson output results ...')
    out_data = dict()
    for provider, provider_gdf in floodings_providers_gdfs.items():
        provider_feature_collection = provider_gdf.to_geo_dict()
        provider_feature_collection['metadata'] = {
            'provider': provider,
            'buildings_count': len(provider_gdf),
            'flooded_buildings_count': len(provider_gdf['is_flooded'])
        }
        provider_feature_collection['crs'] = {
            "type": "name",
            "properties": {
                "name": f"urn:ogc:def:crs:{t_srs.replace(':', '::')}"  # REF: https://gist.github.com/sgillies/1233327 lines 256:271
            }
        }
        print(f"## Provider {provider} feature collection prepared with {len(provider_feature_collection['features'])} features.")
        out_data[provider] = provider_feature_collection
        
    
    # DOC: 7 — Save results to file
    print(f'# Saving results to {out} ...')
    for provider, provider_feature_collection in out_data.items():
        out_provider_fname = f'{out}__{provider}.geojson'
        with open(out_provider_fname, 'w') as f:
            json.dump(provider_feature_collection, f, indent=2)
        print(f"## Results for provider {provider} saved to {out_provider_fname}")
    return out_data



# DOC: Main function to run the flooded buildings analysis from command line

@click.command()
@click.option('--wd', type=click.Path(exists=True), required=True, help='Path to the water depth raster file.')
@click.option('--buildings', type=click.Path(exists=True), default=None, help='Path to the buildings vector file.')
@click.option('--bbox', type=float, nargs=4, default=None, help='Bounding box (minx, miny, maxx, maxy).')
@click.option('--out', type=click.Path(), default=None, help='Output path for the results.')
@click.option('--t_srs', type=str, default=None, help='Target spatial reference system (EPSG code).')
@click.option('--providers', type=str, multiple=True, default=None, help='List of data providers (can be repeated).')
@click.option('--filters', type=str, default=None, help='Filters for providers-features in JSON format.')
def main(wd, buildings, bbox, out, t_srs, providers, filters):
    """
    Main function to run the flooded buildings analysis from command line.
    
    Example:
    
    safer-buildings --wd tests\rimini-wd.tif --providers OVERTURE --filters "{'OVERTURE': [{'subtype':'education', 'class': ['kindergarten','school']}, {'class':'parking'}]}"
    
    in this example, the water depth file is 'tests/rimini-wd.tif', the OVERTURE provider is used, and buildings are filtered is (subtype in ['education'] AND class in ['kindergarten', 'school']) OR (class in ['parking']).
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
    print(f"## Providers: {providers}")
    print(f"## Feature filters: {filters}")
    
    result = compute_flood(
        waterdepth_filename = wd,
        buildings_filename = buildings,
        bbox = tuple(bbox) if bbox else None,
        out = out,
        t_srs = t_srs,
        providers = providers,
        feature_filters = filters
    )
    
    time_end = time.time()
    print(f"# Flooded buildings analysis completed in {time_end - time_start:.2f} seconds.")