import requests
import pandas as pd

# DOC: define costants

_RING_BUFFER_M = 5  # DOC: Default buffer radius in meters for building rings around flooded buildings.

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
    df_layers['provider_name'] = df_layers['id'].apply(lambda service_id: f'RER-REST/{service_id}')
    return df_layers

RegioneEmiliaRomagnaLayers = _list_rest_layers()
_RER_BUILDING_POINTS_BUFFER_M = 20

_PROVIDERS = (
    'OVERTURE',
    * RegioneEmiliaRomagnaLayers.provider_name.to_list(),
)

