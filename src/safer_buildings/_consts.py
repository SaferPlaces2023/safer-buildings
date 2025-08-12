import urllib3
import requests
import xmltodict
import pandas as pd

# DOC: define costants


_RING_BUFFER_M = 5  # DOC: Default buffer radius in meters for building rings around flooded buildings.

_ERROR_CODE_KEY = 'status_code'
_OUTPUT_SUCCESS_CODE = 200


# DOC: List of providers for safer buildings

_PROVIDERS = [ ]


# DOC: Overture provider definition

_OVERTURE_PROVIDER = 'OVERTURE'
_PROVIDERS.append(_OVERTURE_PROVIDER)



# DOC: RER REST provider defintion

_RER_REST_PROVIDER = 'RER-REST'
_RER_BUILDING_POINTS_BUFFER_M = 20
_RER_REST_SERVICE_URL = "https://servizigis.regione.emilia-romagna.it/geoags/rest/services/portale/saferplaces/MapServer"

RegioneEmiliaRomagnaLayers = None

def init_rer_rest_layers():
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
    
    global RegioneEmiliaRomagnaLayers
    RegioneEmiliaRomagnaLayers = pd.DataFrame(data['layers']).sort_values('id').reset_index(drop=True)
    RegioneEmiliaRomagnaLayers['provider_name'] = RegioneEmiliaRomagnaLayers['id'].apply(lambda service_id: f'{_RER_REST_PROVIDER}/{service_id}')
    global _PROVIDERS
    _PROVIDERS.append(_RER_REST_PROVIDER)
    _PROVIDERS.extend(RegioneEmiliaRomagnaLayers['provider_name'].to_list())



# DOC: Venezia WFS provider definition

_VENEZIA_WFS_PROVIDER = 'VENEZIA-WFS'
_VENICE_BUILDING_POINTS_BUFFER_M = 20
_VENEZIA_WFS_SERVICE_URL = "https://webgis2.cittametropolitana.ve.it/lizmap/index.php/lizmap/service?repository=in4safety&project=protezionecivile_ogc&service=WFS&version=1.1.0"

VeneziaLayers = None

def init_venezia_wfs_layers():
    # DOC: Disable SSL warnings for insecure requests (only for VENEZIA WFS that needs verify parameter set to False)    
    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
    params = {
        'request': 'GetCapabilities',
    }
    response = requests.get(_VENEZIA_WFS_SERVICE_URL, params=params, verify=False)
    data = xmltodict.parse(response.content.decode('utf-8'))
    global VeneziaLayers
    VeneziaLayers = pd.DataFrame(data['WFS_Capabilities']['FeatureTypeList']['FeatureType'])
    VeneziaLayers['provider_name'] = VeneziaLayers['Name'].apply(lambda name: f'{_VENEZIA_WFS_PROVIDER}/{name}')
    global _PROVIDERS
    _PROVIDERS.append(_VENEZIA_WFS_PROVIDER)
    _PROVIDERS.extend(VeneziaLayers['provider_name'].to_list())



# DOC: Flood modes definition
class FloodModes:
    """
    Class to define flood modes.
    """
    BUFFER = 'BUFFER'
    IN_AREA = 'IN-AREA'
    ALL = 'ALL'

_FLOOD_MODES = { FloodModes.BUFFER, FloodModes.IN_AREA, FloodModes.ALL }



# DOC: Additional column names created during processing
_COL_FLOOD_ROI = '__flood_roi__'
_COL_FLOOD_AREA = '__flood_area__'
_COL_FLOOD_GEOMETRY = '__flood_geometry__'
_COL_FLOOD_COORDS = '__flood_coords__'
_COL_IS_FLOODED = 'is_flooded'  # DOC: this has not dashes beacuse it will be in the final output as property name

# DOC: Garbage collection of temp-files
_GARBAGE_TEMP_FILES = set()