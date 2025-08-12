import io
import requests
import pandas as pd

from shapely.geometry import box, Point, Polygon, MultiPolygon, LineString, MultiLineString
import geopandas as gpd

import leafmap

from . import _consts, _utils, filesystem

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
        
        if provider == _consts._OVERTURE_PROVIDER:
            provider_buildings = retrieve_overture(bbox)
        elif provider.startswith(_consts._RER_REST_PROVIDER):
            provider_buildings = retrieve_rer_rest(provider, bbox)     
        elif provider.startswith(_consts._VENEZIA_WFS_PROVIDER):
            provider_buildings = retrieve_venezia_wfs(provider, bbox) 
        else:
            raise ValueError(f"Provider '{provider}' is not supported. Available providers are: {_consts._PROVIDERS}.")      
        
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
    service_ids = list(map(int, provider.split('/')[1:])) if provider != 'RER-REST' else _consts.RegioneEmiliaRomagnaLayers[pd.isnull(_consts.RegioneEmiliaRomagnaLayers.subLayerIds)].id.unique().tolist()
    while any([_consts.RegioneEmiliaRomagnaLayers[_consts.RegioneEmiliaRomagnaLayers.id == service_id].iloc[0].subLayerIds is not None for service_id in service_ids]):
        for service_id in service_ids:
            sub_layers = _consts.RegioneEmiliaRomagnaLayers[_consts.RegioneEmiliaRomagnaLayers.id == service_id].iloc[0].subLayerIds
            if sub_layers is not None:
                service_ids.remove(service_id)
                service_ids.extend(sub_layers)
                service_ids = list(set(service_ids))
                
    def rest_service_retrieve(service_id):
        service_class = _consts.RegioneEmiliaRomagnaLayers[_consts.RegioneEmiliaRomagnaLayers.id == service_id].iloc[0]['name']
        url = f"{_consts._RER_REST_SERVICE_URL}/{service_id}/query"
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
        rest_gdf['service_id'] = service_id
        rest_gdf['service_class'] = service_class
        Logger.debug(f"### Retrieved {len(rest_gdf)} features from {service_id} ({service_class}) {_consts._RER_REST_PROVIDER} service.")
        return rest_gdf
    
    def overture_intersection(gdf_re):
        Logger.debug(f"### Overture intersection for {len(gdf_re)} {_consts._RER_REST_PROVIDER} Points.")
        gdf_ot = retrieve_overture(bbox)
        gdf_re['ot_id'] = gdf_re.geometry.apply(lambda geom: gdf_ot[gdf_ot.geometry.contains(geom)].id.values.tolist() if type(geom) is Point else None)
        gdf_re['ot_id'] = gdf_re['ot_id'].apply(lambda ids: ids[0] if ids is not None and len(ids) > 0 else None)
        gdf_re['geometry'] = gdf_re.apply(lambda row: gdf_ot[gdf_ot.id == row.ot_id].iloc[0].geometry if row.ot_id is not None else row.geometry, axis=1)
        Logger.debug(f'### Overture intersection: Taking overture building from {len(gdf_re[gdf_re.ot_id.notnull()])} original {_consts._RER_REST_PROVIDER} Points.')
        return gdf_re

    provider_buildings = pd.concat([rest_service_retrieve(service_id) for service_id in service_ids], ignore_index=True)
    provider_buildings = overture_intersection(provider_buildings)
    provider_buildings = _utils.buffer_points(provider_buildings, buffer_meters=_consts._RER_BUILDING_POINTS_BUFFER_M)
    
    buildings_filename = _utils.temp_filename(ext='gpkg', prefix=f"safer-buildings_{provider.replace('/','-')}")
    provider_buildings.rename(columns={'fid': '_fid'}, inplace=True, errors='ignore')
    provider_buildings.to_file(buildings_filename, driver=filesystem._GPD_DRIVERS('gpkg'), index=False)

    Logger.debug(f"### Retrieved {len(provider_buildings)} features from {provider} service. Saved at {buildings_filename}.")

    return provider_buildings


def retrieve_venezia_wfs(provider, bbox, buffer_points=True):
    service_ids = None
    if provider == _consts._VENEZIA_WFS_PROVIDER:
        service_ids = _consts.VeneziaLayers.Name.unique().tolist()
    elif provider == _consts._VENEZIA_WFS_CRITICAL_SITES_PROVIDER:
        service_ids = _consts.get_venezia_wfs_criticals_layers()
    else:
        service_ids = list(provider.split('/')[1:]) if provider != _consts._VENEZIA_WFS_PROVIDER else _consts.VeneziaLayers.Name.unique().tolist()

    gdf_layers = []
    bounds = bbox.total_bounds
    for service_id in service_ids:
        params={
            "request":"GetFeature",
            "TYPENAME": service_id,
            "outputFormat": "application/json",
            "srsName": _consts.VeneziaLayers[_consts.VeneziaLayers.Name == service_id].iloc[0].DefaultSRS # ???: srsname = _consts.VeneziaLayers[_consts.VeneziaLayers.Name == service_id].iloc[0].DefaultSRS, allow_override=True` but data from server are sent in 4326
        }
        response_data = requests.get(_consts._VENEZIA_WFS_SERVICE_URL, params=params, verify=False)
        geojson_io = io.StringIO(response_data.text)
        wfs_gdf = gpd.read_file(geojson_io).set_crs(crs="EPSG:4326").to_crs(bbox.crs)
        wfs_gdf = wfs_gdf.cx[bounds[0]:bounds[2], bounds[1]:bounds[3]]
        wfs_gdf['service_id'] = service_id
        gdf_layers.append(wfs_gdf)
        Logger.debug(f"### Retrieved {len(wfs_gdf)} features from {service_id} {_consts._VENEZIA_WFS_PROVIDER} service.")

    provider_buildings = pd.concat(gdf_layers, ignore_index=True)
    if buffer_points:
        provider_buildings = _utils.buffer_points(provider_buildings, buffer_meters=_consts._VENEZIA_BUILDING_POINTS_BUFFER_M)

    buildings_filename = _utils.temp_filename(ext='gpkg', prefix=f"safer-buildings_{provider.replace('/','-')}" if len(provider.split('/')) <= 2 else "safer-buildings_venezia-wfs")
    provider_buildings.rename(columns={'fid': '_fid'}, inplace=True, errors='ignore')
    provider_buildings.to_file(buildings_filename, driver=filesystem._GPD_DRIVERS('gpkg'), index=False)

    Logger.debug(f"### Retrieved {len(provider_buildings)} features from {provider} service. Saved at {buildings_filename}.")

    return provider_buildings