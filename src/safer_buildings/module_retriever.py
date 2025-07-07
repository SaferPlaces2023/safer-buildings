import requests
import pandas as pd

from shapely.geometry import box, Point, Polygon, MultiPolygon, LineString, MultiLineString
import geopandas as gpd

import leafmap

from . import _utils
from . import _consts

from .module_log import Logger


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
        bbox = bbox.to_crs(_utils.get_geodataframe_crs(provider_buildings)).total_bounds
        provider_buildings = provider_buildings.cx[bbox[0]:bbox[2], bbox[1]:bbox[3]]
        Logger.debug(f"## Using provided buildings data from {buildings_filename}. Found {len(provider_buildings)} buildings.")
    
    else:
        Logger.debug(f"## Retrieving buildings data from provider: {provider} ...")
            
        buildings_filename = _utils.temp_filename(ext='shp', prefix=f"{provider.replace('/','-')}_buildings")
        
        if provider == 'OVERTURE':
            provider_buildings = retrieve_overture(bbox)
        elif provider.startswith('RER-REST'):
            provider_buildings = retrieve_rer_rest(provider, bbox)      
        else:
            raise ValueError(f"Provider '{provider}' is not supported. Available providers are: {_consts._PROVIDERS}.")

        Logger.debug(f"### Buildings data retrieved from {provider} saved at {buildings_filename}. (Found {len(provider_buildings)} buildings)")        
        
    return provider_buildings


def retrieve_overture(bbox):
    columns = ["id", "geometry", "height", "subtype", "class", "is_underground"]
    provider_buildings = leafmap.get_overture_data(
        overture_type = "building", 
        bbox = bbox.to_crs(epsg=4326).total_bounds.tolist(), 
        columns = columns,
        output = _utils.temp_filename(ext='geojson', prefix='safer-buildings_overture-buildings')
    )
    return provider_buildings


def retrieve_rer_rest(provider, bbox):
    # DOC: Retrieve from all subservices of the requested one.
    service_ids = list(map(int, provider.split('/')[1:]))
    while any([_consts.RegioneEmiliaRomagnaLayers[_consts.RegioneEmiliaRomagnaLayers.id == service_id].iloc[0].subLayerIds is not None for service_id in service_ids]):
        for service_id in service_ids:
            sub_layers = _consts.RegioneEmiliaRomagnaLayers[_consts.RegioneEmiliaRomagnaLayers.id == service_id].iloc[0].subLayerIds
            if sub_layers is not None:
                service_ids.remove(service_id)
                service_ids.extend(sub_layers)
                service_ids = list(set(service_ids))
                
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
        
        geometry_type = _consts.RegioneEmiliaRomagnaLayers[_consts.RegioneEmiliaRomagnaLayers.id == service_id].iloc[0].geometryType
        if geometry_type == 'esriGeometryPoint':
            geometries = [
                Point(f['geometry']['x'], f['geometry']['y']) 
                for f in data['features']
            ]
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
        else:
            raise ValueError(f"Unsupported geometry type for service {service_id}: {_consts.RegioneEmiliaRomagnaLayers[_consts.RegioneEmiliaRomagnaLayers.id == service_id].iloc[0].geometryType}")

        rest_gdf = gpd.GeoDataFrame(features, geometry=geometries, crs=f"EPSG:{data['spatialReference']['wkid']}")
        rest_gdf['rer_id'] = service_id
        rest_gdf['rer_class'] = _consts.RegioneEmiliaRomagnaLayers[_consts.RegioneEmiliaRomagnaLayers.id == service_id].name.iloc[0]
            
        return rest_gdf
    
    def overture_intersection(gdf_re):
        Logger.debug(f"### Overture intersection for {len(gdf_re)} RER-REST Points.")
        gdf_ot = retrieve_overture()
        gdf_re['ot_id'] = gdf_re.geometry.apply(lambda geom: gdf_ot[gdf_ot.geometry.contains(geom)].id.values.tolist() if type(geom) is Point else None)
        gdf_re['ot_id'] = gdf_re['ot_id'].apply(lambda ids: ids[0] if ids is not None and len(ids) > 0 else None)
        gdf_re['geometry'] = gdf_re.apply(lambda row: gdf_ot[gdf_ot.id == row.ot_id].iloc[0].geometry if row.ot_id is not None else row.geometry, axis=1)
        Logger.debug(f'### Overture intersection: Taking overture building from {len(gdf_re[gdf_re.ot_id.notnull()])} original RER-REST Points.')
        return gdf_re

    def buffer_points(gdf, buffer_meters = 20):
        gdf = gdf.to_crs("EPSG:7791")
        gdf['geometry'] = gdf.geometry.apply(lambda g: g.buffer(buffer_meters) if type(g) in [Point, LineString, MultiLineString] else g)
        gdf = gdf.to_crs("EPSG:4326")
        return gdf
    
    provider_buildings = pd.concat([rest_service_retrieve(service_id) for service_id in service_ids], ignore_index=True)
    provider_buildings = overture_intersection(provider_buildings)
    provider_buildings = buffer_points(provider_buildings, buffer_meters=_consts._RER_BUILDING_POINTS_BUFFER_M)
    
    buildings_filename = _utils.temp_filename(ext='shp', prefix=f"safer-buildings_{provider.replace('/','-')}")
    provider_buildings.to_file(buildings_filename, driver='ESRI Shapefile', index=False)

    return provider_buildings