import requests
import xmltodict
import pandas as pd

# DOC: define costants


_RING_BUFFER_M = 5  # DOC: Default buffer radius in meters for building rings around flooded buildings.


_RER_REST_SERVICE_URL = "https://servizigis.regione.emilia-romagna.it/geoags/rest/services/portale/saferplaces/MapServer"

def _rer_list_rest_layers():
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
    response = requests.get(_RER_REST_SERVICE_URL, params=params, headers=headers)
    data = response.json()
    df_layers = pd.DataFrame(data['layers']).sort_values('id').reset_index(drop=True)
    df_layers['provider_name'] = df_layers['id'].apply(lambda service_id: f'RER-REST/{service_id}')
    return df_layers

RegioneEmiliaRomagnaLayers = _rer_list_rest_layers()
_RER_BUILDING_POINTS_BUFFER_M = 20


_VENEZIA_WFS_SERVICE_URL = "https://webgis2.cittametropolitana.ve.it/lizmap/index.php/lizmap/service?repository=in4safety&project=gecosistema_ogc&service=WFS&version=1.1.0"

def _venezia_list_wfs_layers():
    params = {
        'request': 'GetCapabilities',
    }
    response = requests.get(_VENEZIA_WFS_SERVICE_URL, params=params, verify=False)
    data = xmltodict.parse(response.content.decode('utf-8'))
    df_layers = pd.DataFrame(data['WFS_Capabilities']['FeatureTypeList']['FeatureType'])
    df_layers['provider_name'] = df_layers['Name'].apply(lambda name: f'VENEZIA-WFS/{name}')
    return df_layers

VeneziaLayers = _venezia_list_wfs_layers()
_VENICE_BUILDING_POINTS_BUFFER_M = 20



_PROVIDERS = (
    'OVERTURE',

    'RER-REST',
    * RegioneEmiliaRomagnaLayers.provider_name.to_list(),

    'VENEZIA-WFS',
    * VeneziaLayers.provider_name.to_list(),
)

