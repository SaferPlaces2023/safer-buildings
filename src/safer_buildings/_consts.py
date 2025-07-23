import requests
import xmltodict
import pandas as pd

# DOC: define costants


_RING_BUFFER_M = 5  # DOC: Default buffer radius in meters for building rings around flooded buildings.



# DOC: Overture provider definition

_OVERTURE_PROVIDER = 'OVERTURE'



# DOC: RER REST provider defintion

_RER_REST_PROVIDER = 'RER-REST'
_RER_BUILDING_POINTS_BUFFER_M = 20
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
    df_layers['provider_name'] = df_layers['id'].apply(lambda service_id: f'{_RER_REST_PROVIDER}/{service_id}')
    return df_layers

RegioneEmiliaRomagnaLayers = _rer_list_rest_layers()



# DOC: Venezia WFS provider definition

_VENEZIA_WFS_PROVIDER = 'VENEZIA-WFS'
_VENICE_BUILDING_POINTS_BUFFER_M = 20
_VENEZIA_WFS_SERVICE_URL = "https://webgis2.cittametropolitana.ve.it/lizmap/index.php/lizmap/service?repository=in4safety&project=gecosistema_ogc&service=WFS&version=1.1.0"

def _venezia_list_wfs_layers():
    params = {
        'request': 'GetCapabilities',
    }
    response = requests.get(_VENEZIA_WFS_SERVICE_URL, params=params, verify=False)
    data = xmltodict.parse(response.content.decode('utf-8'))
    df_layers = pd.DataFrame(data['WFS_Capabilities']['FeatureTypeList']['FeatureType'])
    df_layers['provider_name'] = df_layers['Name'].apply(lambda name: f'{_VENEZIA_WFS_PROVIDER}/{name}')
    return df_layers

VeneziaLayers = _venezia_list_wfs_layers()



# DOC: List of providers for safer buildings

_PROVIDERS = (
    _OVERTURE_PROVIDER,

    _RER_REST_PROVIDER,
    * RegioneEmiliaRomagnaLayers.provider_name.to_list(),

    _VENEZIA_WFS_PROVIDER,
    * VeneziaLayers.provider_name.to_list(),
)

