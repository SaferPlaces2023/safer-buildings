import urllib3
import requests
import xmltodict
import pandas as pd

# DOC: define costants


_RING_BUFFER_M = 5  # DOC: Default buffer radius in meters for building rings around flooded buildings.

_ERROR_CODE_KEY = 'status_code'
_OUTPUT_SUCCESS_CODE = 200


# DOC: Projected EPSG for spatial operations
_EPSG_UTMxx = None
def init_epsg_utmxx(epsg_code):
    """
    Set the EPSG code for UTM projection.
    """
    global _EPSG_UTMxx
    _EPSG_UTMxx = epsg_code


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
_VENEZIA_WFS_CRITICAL_SITES_PROVIDER = 'VENEZIA-WFS-CRITICAL-SITES'
_VENEZIA_BUILDING_POINTS_BUFFER_M = 20
_VENEZIA_WFS_SERVICE_URL = "https://webgis2.cittametropolitana.ve.it/lizmap/index.php/lizmap/service?repository=in4safety&project=protezionecivile_ogc&service=WFS&version=1.1.0"

_VENEZIA_WFS_NAME2LABEL_MAP = labels = {
    'c0107057_grafostrade': 'Grafo strade',
    'mv_risorse_p0109033_materiali': 'Materiali',
    'mv_risorse_p0109043_natanti': 'Natanti',
    'mv_risorse_p0109053_mezziaerei': 'Mezzi aerei',
    'mv_risorse_p0109063_mezziautomotoveicoli': 'Mezzi auto-motoveicoli',
    'mv_risorse_p0109073_macchineoperatrici': 'Macchine operatrici',
    'mv_risorse_p0109083_carrellielevatori': 'Carrelli elevatori',
    'mv_risorse_p0109093_rimorchi': 'Rimorchi',
    'mv_risorse_p0109103_pompe': 'Pompe',
    'mv_risorse_p0109113_generatori': 'Generatori',
    'mv_risorse_p0109123_fari': 'Fari',
    'mv_risorse_p0109133_moduli': 'Moduli',
    'mv_risorse_p0109143_cucine': 'Cucine',
    'mv_risorse_p0109153_cisterne': 'Cisterne',
    'mv_risorse_p0109163_tende': 'Tende',
    'mv_risorse_p0109183_altreattrezzature': 'Altre attrezzature',
    'v_pc_p0102011_areeattesa': 'Aree di attesa',
    'v_pc_p0102021_areericovero': 'Aree di ricovero',
    'v_pc_p0102031_areeammassamento': 'Aree di ammassamento',
    'v_pc_p0103011_allertamento': 'Aree allertamento',
    'v_pc_p0103013_allertamento': 'Metodi allertamento',
    'v_pc_p0104011_ufficipc': 'Uffici di Protezione Civile',
    'v_pc_p0104021_magazzinipc': 'Magazzini di Protezione Civile',
    'v_pc_p0104031_prontosoccorso': 'Pronto soccorso',
    'v_pc_p0104041_sediamministrative': 'Sedi amministrative',
    'v_pc_p0104051_strutturestrategichespecifiche': 'Strutture strategiche specifiche',
    'v_pc_p0105011_stazioniferroviarie': 'Stazioni ferroviarie',
    'v_pc_p0105021_busmetro': 'Bus e metro',
    'v_pc_p0105031_aeroporti': 'Aeroporti',
    'v_pc_p0105041_porti': 'Porti',
    'v_pc_p0105052_stradestrategiche': 'Strade strategiche',
    'v_pc_p0105062_pontistrategici': 'Ponti strategici',
    'v_pc_p0105072_galleriestrategiche': 'Gallerie strategiche',
    'v_pc_p0105081_operepresa': 'Opere di presa',
    'v_pc_p0105092_acquedotti': 'Acquedotti',
    'v_pc_p0105101_trasformazioneee': 'Trasformazione energia elettrica',
    'v_pc_p0105112_elettrodotti': 'Elettrodotti',
    'v_pc_p0105121_news': 'Notizie',
    'v_pc_p0105131_distributoricarburante': 'Distributori di carburante',
    'v_pc_p0106011_scuole': 'Scuole',
    'v_pc_p0106021_stadi': 'Stadi',
    'v_pc_p0106031_strutturesanitarie': 'Strutture sanitarie',
    'v_pc_p0106041_casecura': 'Case di cura',
    'v_pc_p0106051_edificirilevantigenerici': 'Edifici rilevanti (generici)',
    'v_pc_p0106061_centricommerciali': 'Centri commerciali',
    'v_pc_p0106071_museibiblioteche': 'Musei e biblioteche',
    'v_pc_p0106081_salespettacoli': 'Sale spettacoli',
    'v_pc_p0106091_banchept': 'Banche e uffici postali',
    'v_pc_p0106101_industrie': 'Industrie',
    'v_pc_p0106111_alberghi': 'Alberghi',
    'v_pc_p0106121_localiculto': 'Luoghi di culto',
    'v_pc_p0107012_collegamentiviari': 'Collegamenti viari',
    'v_pc_p0107022_ponti': 'Ponti',
    'v_pc_p0107032_gallerie': 'Gallerie',
    'v_pc_p0107042_dighe': 'Dighe',
    'v_pc_p0108011_prefetture': 'Prefetture',
    'v_pc_p0108021_depositofarmaci': 'Depositi farmaci',
    'v_pc_p0108031_depositoalimenti': 'Depositi alimenti',
    'v_pc_p0108041_allevamenti': 'Allevamenti',
    'v_pc_p0108051_carceri': 'Carceri',
    'v_pc_p0108061_strutturemilitari': 'Strutture militari',
    'v_pc_p0108072_nodisensibili': 'Nodi sensibili',
    'v_pc_p0108081_depuratori': 'Depuratori',
    'v_pc_p0108091_discariche': 'Discariche',
    'v_pc_p0108103_cancelli': 'Cancelli',
    'v_pc_p0108111_cimiteri': 'Cimiteri',
    'v_pc_p0109011_orgvolontariato': 'Organizzazioni di volontariato',
    'v_pc_p0201011_sisma': 'Sisma',
    'v_pc_p0201021_blackout': 'Blackout',
    'v_pc_p0201032_neve': 'Neve',
    'v_pc_p0201042_incidentistradali': 'Incidenti stradali',
    'v_pc_p0201043_accessipma': 'Accessi PMA',
    'v_pc_p0201051_incidentirilevanti': 'Incidenti rilevanti',
    'v_pc_p0201061_zoneimpatto': 'Zone di impatto',
    'v_pc_p0201072_trasportopericolose': 'Trasporti di merci pericolose',
    'v_pc_p0201081_allagamenti': 'Allagamenti',
    'v_pc_p0201091_frane': 'Frane',
    'v_pc_p0201101_mareggiate': 'Mareggiate',
    'v_pc_p0201111_valanghe': 'Valanghe',
    'v_pc_p0201121_crollodighe': 'Crollo dighe',
    'v_pc_p0201131_idropotabile': 'Idropotabile',
    'v_pc_p0202013_idranti': 'Idranti',
    'v_pc_p0202022_stradeforestali': 'Strade forestali',
    'v_pc_p0202032_ostacolivolo': 'Ostacoli al volo'
}


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
    VeneziaLayers['service_label'] = VeneziaLayers['Name'].apply(lambda name: _VENEZIA_WFS_NAME2LABEL_MAP.get(name, name))
    global _PROVIDERS
    _PROVIDERS.append(_VENEZIA_WFS_PROVIDER)
    _PROVIDERS.extend(VeneziaLayers['provider_name'].to_list())

def get_venezia_wfs_criticals_layers():
    if VeneziaLayers is None:
        init_venezia_wfs_layers()
    critical_layers_titles = [
        title for title in VeneziaLayers.Title 
        if title.startswith('p0106') or title=='p0104031_prontosoccorso' or title=='c0107057_grafostrade'
    ]
    critical_layers_names = VeneziaLayers[VeneziaLayers.Title.isin(critical_layers_titles)].Name.unique().tolist()
    return critical_layers_names



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