import os
import json
import datetime
from tqdm import tqdm

import pandas as pd
import geopandas as gpd
from shapely import distance
from shapely.geometry import box, MultiPoint, MultiPolygon
from shapely.strtree import STRtree

from . import _consts, _utils
from . import module_retriever
from .module_log import Logger, is_debug_mode



class AdditionalOperation():

    def __init__(self, name: str):
        self.name = name

    def __str__(self):
        return f"Callable: {self.name}"
        



# DOC: This is for VENEZIA_WFS_PROVIDER (Find the nearest pump to a flooded building)
class NearbyPumps(AdditionalOperation):
    
    """
    This operation retrieves the nearby pumps for flooded buildings based on water depth areas. 
    It finds pumps within a specified distance from the water depth areas.
    The operation uses the `mv_risorse_p0109103_pompe` layer from the VENEZIA WFS provider to find pumps.
    The avaliable pumps for each flooded area are then associated with the related flooded buildings, and the distance to each pump from them is calculated
    """

    name = 'nearby_pumps'
    
    description = {
        'en': """This operation retrieves the nearby pumps for flooded buildings based on water depth areas. 
        It finds pumps within a specified distance from the water depth areas.
        The operation uses the `mv_risorse_p0109103_pompe` layer from the VENEZIA WFS provider to find pumps.
        The avaliable pumps for each flooded area are then associated with the related flooded buildings, and the distance to each pump from them is calculated.""",
        
        'it': """Questa operazione recupera le pompe vicine per gli edifici allagati basandosi sulle aree di profondità dell'acqua.
        Trova le pompe entro una distanza specificata dalle aree di profondità dell'acqua.
        L'operazione utilizza il layer `mv_risorse_p0109103_pompe` dal provider VENEZIA WFS per trovare le pompe.
        Le pompe disponibili per ciascuna area allagata vengono quindi associate agli edifici allagati correlati, e viene calcolata la distanza da ciascuna pompa.""",
    }
    
    args = [
        # DOC: The maximum distance to consider a pump as nearby (in meters). Default is 1000.0 meters.
        'max_distance'          
    ]

    _pumps_layer_id = 'mv_risorse_p0109103_pompe'
    _pumps_identifier = 'gid'
    _nearby_pumps_basic_attributes = ['id', '_fid', 'gid', 'modello', 'indirizzo']

    def __init__(self, max_distance: float = 1000.0):
        super().__init__(name=self.name)
        self._configure(max_distance=max_distance)

    def _configure(self, max_distance: float = 100.0):
        self.max_distance = float(max_distance) if type(max_distance) is str else max_distance

    def __call__(self, **kwargs):

        # DOC: Extract kwargs values ----------------------------------------------------------------------------------
        gdf_buildings = kwargs['gdf_buildings'].reset_index(drop=True)
        gdf_wd = kwargs['gdf_water_depth']
        bbox = kwargs['bbox']        
        t_srs = kwargs['t_srs']
        
        # DOC: Convert CRS to Projected CRS (EPSG:UTMxx) due to calculating distances
        gdf_buildings['geometry_prj'] = gdf_buildings.to_crs(_consts._EPSG_UTMxx).geometry

        # DOC: Retrieve the pumps from the WFS provider ---------------------------------------------------------------
        gdf_pumps = module_retriever.retrieve_venezia_wfs(
            provider = f'{_consts._VENEZIA_WFS_PROVIDER}/{self._pumps_layer_id}',
            bbox = bbox,
            buffer_points = False 
        ).reset_index(drop=True)
        # DOC: Prepare a projected-crs geometry column to use in spatial operations
        gdf_pumps['geometry_prj'] = gdf_pumps.to_crs(_consts._EPSG_UTMxx).geometry
        gdf_pumps['geometry'] = gdf_pumps['geometry_prj'].centroid
        # DOC: Convert the pumps GeoDataFrame to the target CRS
        gdf_pumps = gdf_pumps.to_crs(t_srs)
        gdf_pumps['longitude'] = gdf_pumps.geometry.x
        gdf_pumps['latitude'] = gdf_pumps.geometry.y
        self._nearby_pumps_basic_attributes.extend(['longitude', 'latitude'])

        # DOC: Create a spatial index for the water depth polygons ----------------------------------------------------
        wd_prj_tree = STRtree(gdf_wd.to_crs(_consts._EPSG_UTMxx).geometry.values)

        # DOC: Find pumps whithin a distance from the water depth polygons --------------------------------------------
        pumps_in_wd_query = wd_prj_tree.query(gdf_pumps['geometry_prj'], predicate='dwithin', distance=self.max_distance)
        pumps_selection = gdf_pumps.iloc[pumps_in_wd_query[0,:]]   # DOC: Extract from Target CRS pumps df, using the query result, so we have the pumps in the target CRS
        if pumps_selection.empty:
            Logger.debug("## No pumps found within the specified distance from water depth areas.")
            return gdf_buildings, { 'type': 'FeatureCollection', 'features': [] }
        # DOC: Prepare a DataFrame with the pumps and their water depth area index
        pumps_selection['index'] = pumps_selection.index
        pumps_selection = pumps_selection.drop_duplicates(subset='index').drop(columns='index')
        gdf_pumps = pd.merge(
            pumps_selection, pd.DataFrame(pumps_in_wd_query.T, columns=['pump_idx', 'wd_idx']),
            left_index=True, right_on='pump_idx', how='left'
        ).drop(columns='pump_idx')
        gdf_pumps.set_index('wd_idx', inplace=True)
        Logger.debug(f"## Found {len(gdf_pumps[self._pumps_identifier].unique())} nearby pumps for the flooded areas.")

        # DOC: Find wd area for each building -------------------------------------------------------------------------
        buildings_in_wd_query = wd_prj_tree.query(gdf_buildings['geometry_prj'], predicate='intersects', distance=0)
        bld_wd_idxs = pd.DataFrame(buildings_in_wd_query.T, columns=['building_idx', 'wd_idxs']).groupby('building_idx').agg(list).reset_index()

        # DOC: Join the buildings with related water depth areas indexes ----------------------------------------------
        gdf_buildings = pd.merge(
            gdf_buildings, bld_wd_idxs,
            left_index=True, right_on='building_idx', how='left'
        ).drop(columns='building_idx')
        # DOC: Explode the buildings GeoDataFrame to have one row per building (duplicates) and water depth area (become single)
        gdf_buildings_pumps = gdf_buildings.explode('wd_idxs')
        # DOC: Join the buildings with pumps that share the same water depth area
        gdf_pumps.rename(columns={col: f'{col}_pump' for col in gdf_pumps.columns}, inplace=True)
        gdf_buildings_pumps = pd.merge(
            gdf_buildings_pumps, gdf_pumps[[f'{attr}_pump' for attr in self._nearby_pumps_basic_attributes] + ['geometry_prj_pump']],
            left_on='wd_idxs', right_index=True, how='left'
        ).drop(columns='wd_idxs')
        gdf_buildings_pumps = gdf_buildings_pumps[gdf_buildings_pumps['geometry_prj_pump'].notna()]
        
        # DOC: Calculate the distance between each building and its related pumps (then filter by max_distance)
        gdf_buildings_pumps['distance_pump'] = distance(gdf_buildings_pumps['geometry_prj'].values, gdf_buildings_pumps['geometry_prj_pump'].values)
        gdf_buildings_pumps = gdf_buildings_pumps[gdf_buildings_pumps['distance_pump'] <= self.max_distance]
        
        # DOC: Foreach building, collect the nearby pumps (only basic attributes) in a list
        gdf_buildings_pumps = gdf_buildings_pumps[[col for col in gdf_buildings_pumps.columns if col.endswith('_pump')]].drop(columns=['geometry_prj_pump'])
        gdf_buildings_pumps.rename(columns={col: col.replace('_pump', '') for col in gdf_buildings_pumps.columns}, inplace=True)
        gdf_buildings_pumps = gdf_buildings_pumps.groupby(gdf_buildings_pumps.index).agg(list)
        buildings_nearby_pumps = gdf_buildings_pumps.apply(lambda row: pd.DataFrame(row.to_dict()).to_dict(orient='records'), axis=1).to_list()
        # DOC: Add the nearby pumps to the buildings GeoDataFrame       
        gdf_buildings['nearby_pumps'] = [list() for _ in range(len(gdf_buildings))]
        gdf_buildings.loc[gdf_buildings_pumps.index, 'nearby_pumps'] = pd.Series(buildings_nearby_pumps, index=gdf_buildings_pumps.index)
        gdf_buildings.drop(columns=['geometry_prj'], inplace=True)

        # DOC: Prepare the pumps GeoDataFrame to be returned ----------------------------------------------------------
        gdf_pumps.rename(columns={col: col.replace('_pump', '') for col in gdf_pumps.columns}, inplace=True)
        gdf_pumps.drop_duplicates(subset=self._pumps_identifier, inplace=True)
        gdf_pumps = gdf_pumps.drop(columns=['geometry_prj']).reset_index(drop=True)
        gdf_pumps = _utils.safe_json_df(gdf_pumps)
        nearby_pumps_collection = gdf_pumps.to_geo_dict()
        nearby_pumps_collection = _utils.set_crs_feature_collection(nearby_pumps_collection, t_srs)
        Logger.debug(f"## Buildings with nearby pumps: {gdf_buildings[gdf_buildings[_consts._COL_IS_FLOODED]]['nearby_pumps'].apply(lambda npumps: len(npumps)>0).sum()} out of {len(gdf_buildings[gdf_buildings[_consts._COL_IS_FLOODED]])} total flooded buildings.")

        # DOC: Return the updated GeoDataFrame with the nearest pumps and the nearby pumps collection -----------------
        return gdf_buildings, nearby_pumps_collection
            



# DOC: This is for VENEZIA_WFS_PROVIDER (Alert method to use)
class AlertMethod(AdditionalOperation):

    """
    This operation retrieves the alert methods for flooded buildings based on water depth areas.
    It finds alert methods within a specified buffer around the water depth areas.
    The operation uses the `v_pc_p0103011_allertamento` layer from the VENEZIA WFS provider to find alert areas,
    and the `v_pc_p0103013_allertamento` layer to find alert methods.
    The avaliable alert methods for each flooded area are then associated with the related flooded buildings.
    """

    name = 'alert_method'

    description = {
        'en': """This operation retrieves the alert methods for flooded buildings based on water depth areas.
        It finds alert methods within a specified buffer around the water depth areas.
        The operation uses the `v_pc_p0103011_allertamento` layer from the VENEZIA WFS provider to find alert areas,
        and the `v_pc_p0103013_allertamento` layer to find alert methods.
        The avaliable alert methods for each flooded area are then associated with the related flooded buildings.""",

        'it': """Questa operazione recupera i metodi di allerta per gli edifici allagati basandosi sulle aree di profondità dell'acqua.
        Trova i metodi di allerta entro un buffer specificato attorno alle aree di profondità dell'acqua.
        L'operazione utilizza il layer `v_pc_p0103011_allertamento` dal provider VENEZIA WFS per trovare le aree di allerta,
        e il layer `v_pc_p0103013_allertamento` per trovare i metodi di allerta.
        I metodi di allerta disponibili per ciascuna area allagata vengono quindi associati agli edifici allagati correlati.""",
    }

    args = [
        # DOC: The buffer to use around water depth areas (in meters). Default is 100.0 meters.
        'wd_buffer'
    ]

    _alert_area_layer_id = 'v_pc_p0103011_allertamento'
    _alert_area_identifier = '_fid'
    _alert_method_layer_id = 'v_pc_p0103013_allertamento'
    _alert_method_identifier = '_fid'
    _alert_method_basic_attributes = ['id', '_fid', 'denom', 'indirizzo', 'strumento_t']

    def __init__(self, wd_buffer: float = 100.0):
        super().__init__(name=self.name)
        self._configure(wd_buffer=wd_buffer)

    def _configure(self, wd_buffer: float = 100.0):
        self.wd_buffer = float(wd_buffer) if type(wd_buffer) is str else wd_buffer
       

    def __call__(self, **kwargs):

        # DOC: Extract kwargs values ----------------------------------------------------------------------------------
        gdf_buildings = kwargs['gdf_buildings'].reset_index(drop=True)
        gdf_wd = kwargs['gdf_water_depth']
        bbox = kwargs['bbox']
        t_srs = kwargs['t_srs']

        # DOC: Convert CRS to Projected CRS (EPSG:UTMxx) due to compute spatial operations
        gdf_buildings['geometry_prj'] = gdf_buildings.to_crs(_consts._EPSG_UTMxx).geometry

        # DOC: Retrieve the alert areas -------------------------------------------------------------------------------
        gdf_alert_area = module_retriever.retrieve_venezia_wfs(
            provider = f'{_consts._VENEZIA_WFS_PROVIDER}/{self._alert_area_layer_id}',
            bbox = bbox,
            buffer_points=False 
        ).reset_index(drop=True)
        gdf_alert_area['geometry_prj'] = gdf_alert_area.to_crs(_consts._EPSG_UTMxx).geometry
        # DOC: Retrieve the alert methods (only related to flooding event) from the WFS provider
        gdf_alert_method = module_retriever.retrieve_venezia_wfs(
            provider = f'{_consts._VENEZIA_WFS_PROVIDER}/{self._alert_method_layer_id}',
            bbox = bbox,
            buffer_points = False
        )
        gdf_alert_method = gdf_alert_method[(gdf_alert_method['allagament_t']=='Sì') | (gdf_alert_method['r_altro']=='Acqua Alta')].reset_index(drop=True)
        gdf_alert_method['geometry_prj'] = gdf_alert_method.to_crs(_consts._EPSG_UTMxx).geometry
        # DOC: Convert the alert method GeoDataFrame to the target CRS
        gdf_alert_method['geometry'] = gpd.GeoSeries(gdf_alert_method['geometry_prj'].centroid, crs=_consts._EPSG_UTMxx).to_crs(gdf_alert_method.crs)
        gdf_alert_method = gdf_alert_method.to_crs(t_srs)
        gdf_alert_method['longitude'] = gdf_alert_method.geometry.x
        gdf_alert_method['latitude'] = gdf_alert_method.geometry.y
        self._alert_method_basic_attributes.extend(['longitude', 'latitude'])

        # DOC: Filter out alert area that not contains alert methods
        gdf_alert_area = gdf_alert_area[gdf_alert_area['geometry_prj'].intersects(gdf_alert_method['geometry_prj'].union_all())].reset_index(drop=True)
        
        # DOC: Create a spatial index for the water depth polygons ----------------------------------------------------
        wd_tree = STRtree(gdf_wd.to_crs(_consts._EPSG_UTMxx).buffer(max(self.wd_buffer, 0)).geometry.values)

        # DOC: Alert area in water depth areas ------------------------------------------------------------------------
        alert_area_in_wd_query = wd_tree.query(gdf_alert_area['geometry_prj'].values, predicate='intersects', distance=0)
        alert_area_selection = gdf_alert_area.iloc[alert_area_in_wd_query[0,:]]
        alert_area_selection['index'] = alert_area_selection.index
        alert_area_selection.drop_duplicates(subset='index', inplace=True)
        alert_area_selection.drop(columns='index', inplace=True)
        Logger.debug(f"## Found {len(alert_area_selection)} alert areas in the water depth areas.")

        # DOC: Alert method in alert area -----------------------------------------------------------------------------
        alert_area_selection_tree = STRtree(alert_area_selection['geometry_prj'].values)
        alert_method_query = alert_area_selection_tree.query(gdf_alert_method['geometry_prj'].values, predicate='within')
        alert_method_selection = gdf_alert_method.iloc[alert_method_query[0,:]]
        alert_method_selection['index'] = alert_method_selection.index
        alert_method_selection.drop_duplicates(subset='index', inplace=True)
        alert_method_selection.drop(columns='index', inplace=True)
        Logger.debug(f"## Found {len(alert_method_selection)} alert methods in the alert areas.")
        # DOC: Join the alert method with the alert area indexes
        gdf_alert_method = pd.merge(
            alert_method_selection, pd.DataFrame({'alert_method_idx': alert_method_query[0,:], 'alert_area_idx': alert_method_query[1,:]}),
            left_index=True, right_on='alert_method_idx', how='left'
        )
        gdf_alert_method.set_index('alert_area_idx', inplace=True)

        # DOC: Prepare and alert method column in the buildings GeoDataFrame ------------------------------------------
        gdf_buildings['alert_method'] = [list() for _ in range(len(gdf_buildings))]
        # DOC: Foreach building, find the alert areas it is in
        buildings_alert_query = alert_area_selection_tree.query(gdf_buildings['geometry_prj'].values, predicate='intersects')
        # DOC: Aggregate alert area indexes for each building
        bld_aa_idxs = pd.DataFrame(buildings_alert_query.T, columns=['building_idx', 'alert_area_idx']).groupby('building_idx').agg(list).reset_index()
        bld_aa_idxs.set_index('building_idx', inplace=True)
        
        # DOC: Join the buildings with related alert area indexes -----------------------------------------------------
        gdf_buildings = pd.merge(
            gdf_buildings, bld_aa_idxs,
            left_index=True, right_index=True, how='left'
        )
        # DOC: Explode the buildings GeoDataFrame to have one row per building (duplicates) and alert area (become single)
        gdf_buildings_alert = gdf_buildings.explode('alert_area_idx')
        # DOC: Join the buildings with alert methods that share the same alert area
        gdf_alert_method.rename(columns={col: f'{col}_alert_method' for col in gdf_alert_method.columns}, inplace=True)
        gdf_buildings_alert = pd.merge(
            gdf_buildings_alert, gdf_alert_method[[f'{attr}_alert_method' for attr in self._alert_method_basic_attributes] + ['geometry_prj_alert_method']],
            left_on='alert_area_idx', right_index=True, how='left'
        ).drop(columns='alert_area_idx')
        gdf_buildings_alert = gdf_buildings_alert[gdf_buildings_alert['geometry_prj_alert_method'].notna()]

        # DOC: Foreach building, collect the alert methods (only basic attributes) in a list
        gdf_buildings_alert = gdf_buildings_alert[[col for col in gdf_buildings_alert.columns if col.endswith('_alert_method')]].drop(columns=['geometry_prj_alert_method'])
        gdf_buildings_alert.rename(columns={col: col.replace('_alert_method', '') for col in gdf_buildings_alert.columns}, inplace=True)
        gdf_buildings_alert = gdf_buildings_alert.groupby(gdf_buildings_alert.index).agg(list)
        buildings_alert_methods = gdf_buildings_alert.apply(lambda row: pd.DataFrame(row.to_dict()).to_dict(orient='records'), axis=1).to_list()
        # DOC: Add the alert methods to the buildings GeoDataFrame
        gdf_buildings['alert_method'] = [list() for _ in range(len(gdf_buildings))]
        gdf_buildings.loc[gdf_buildings_alert.index, 'alert_method'] = pd.Series(buildings_alert_methods, index=gdf_buildings_alert.index)
        gdf_buildings.drop(columns=['geometry_prj'], inplace=True)

        # DOC: Prepare the alert method GeoDataFrame to be returned -------------------------------------------------
        gdf_alert_method.rename(columns={col: col.replace('_alert_method', '') for col in gdf_alert_method.columns}, inplace=True)
        gdf_alert_method.drop_duplicates(subset=self._alert_method_identifier, inplace=True)
        gdf_alert_method = gdf_alert_method.drop(columns=['geometry_prj', 'alert_method_idx']).reset_index(drop=True)
        gdf_alert_method = _utils.safe_json_df(gdf_alert_method)
        alert_method_collection = gdf_alert_method.to_geo_dict()
        alert_method_collection = _utils.set_crs_feature_collection(alert_method_collection, t_srs)
        Logger.debug(f"## Buildings with alert methods: {gdf_buildings[gdf_buildings[_consts._COL_IS_FLOODED]]['alert_method'].apply(lambda am: len(am)>0).sum()} out of {len(gdf_buildings[gdf_buildings[_consts._COL_IS_FLOODED]])} total flooded buildings.")

        # DOC: Return the updated GeoDataFrame with the alert method and the alert method collection ----------------
        return gdf_buildings, alert_method_collection
            



# DOC: This is for VENEZIA_WFS_PROVIDER (Alert method to use)
class GatesGuard(AdditionalOperation):
    """
    This operation retrieves the gates guard for flooded buildings based on flooded streets.
    # TODO: continue ...
    """

    name = 'gates_guard'

    description = {
        'en': """ to be written ... """,

        'it': """ da scrivere ...""",
    }

    args = [
        # DOC: The buffer to use around water depth areas (in meters). Default is 100.0 meters.
        'street_layer_id',
        'street_buffer',
        'flood_area_thresh',
        'gates_buffer',
        'max_distance',
    ]

    _street_graph_layer_id = f'{_consts._VENEZIA_WFS_PROVIDER}/c0107057_grafostrade'
    _gates_layer_id = f'{_consts._VENEZIA_WFS_PROVIDER}/v_pc_p0108103_cancelli'
    
    _street_graph_fid = '_fid'
    _gate_fid = '_fid'
    _gate_basic_attributes = ['id', '_fid', 'denom', 'indirizzo', 'addetti_t']

    _default_street_layer_id = f'{_consts._VENEZIA_WFS_PROVIDER}/v_pc_p0105052_stradestrategiche'  # DOC: other can be `c0107057_grafostrade` but it is huge
    _default_street_buffer = 3.0
    _default_flood_area_thresh = 25.0
    _default_gates_buffer = 3.0
    _default_max_distance = 1000.0  # DOC: max distance to search for gates from flooded street

    def __init__(self, street_layer_id: str = _default_street_layer_id, street_buffer: float = _default_street_buffer, flood_area_thresh: float = _default_flood_area_thresh, gates_buffer: float = _default_gates_buffer, max_distance: float = 5000.0):
        super().__init__(name=self.name)
        self._configure(street_layer_id=street_layer_id, street_buffer=street_buffer, flood_area_thresh=flood_area_thresh, gates_buffer=gates_buffer, max_distance=max_distance)

    def _configure(self, street_layer_id: str = _default_street_layer_id, street_buffer: float = _default_street_buffer, flood_area_thresh: float = _default_flood_area_thresh, gates_buffer: float = _default_gates_buffer, max_distance: float = _default_max_distance):
        if _consts.VeneziaLayers is None:
            _consts.init_venezia_wfs_layers()
        if street_layer_id not in _consts._PROVIDERS:
            raise ValueError(f"Invalid street layer: {street_layer_id}. Available layers: {_consts._PROVIDERS}")
        self.street_layer_id = street_layer_id
        self.street_buffer = float(street_buffer) if type(street_buffer) is str else street_buffer if type(street_buffer) is not None else self._default_street_buffer
        self.flood_area_thresh = float(flood_area_thresh) if type(flood_area_thresh) is str else flood_area_thresh if type(flood_area_thresh) is not None else self._default_flood_area_thresh
        self.gates_buffer = float(gates_buffer) if type(gates_buffer) is str else gates_buffer if type(gates_buffer) is not None else self._default_gates_buffer
        self.max_distance = float(max_distance) if type(max_distance) is str else max_distance if type(max_distance) is not None else self._default_max_distance


    def bfs_streets_gates(self, gdf_streets_graph, street_graph_fid, gdf_gate, max_distance=_default_max_distance):
        """
        BFS to search for gates in the street graph starting from a flooded street.
        """
        streets_graph_fids = [ street_graph_fid ]
        street_gates_df = pd.DataFrame()
        i_start = 0

        street_distances = dict()
        street_distances[street_graph_fid] = 0       

        while len(streets_graph_fids) - i_start > 0:
            s_fid = streets_graph_fids[i_start]
            i_start += 1

            street_graph = gdf_streets_graph[gdf_streets_graph[self._street_graph_fid] == s_fid].iloc[0]
            
            # DOC: If current street has gates the add them to the gate_ids set
            if len(street_graph.gate_idx) > 0:
                found_gates = gdf_gate.loc[street_graph.gate_idx][self._gate_basic_attributes]
                found_gates['distance'] = street_distances[s_fid]
                street_gates_df = pd.concat([street_gates_df,  found_gates], ignore_index=True).reset_index(drop=True)
                del street_distances[s_fid]
                continue
            
            # DOC: Otherwise, find adjacent streets and add them to the queue if not already present
            for adj_fid in filter(lambda adj_fid: adj_fid not in streets_graph_fids, street_graph[f'adjacent_{self._street_graph_fid}']):
                street_graph_adj = gdf_streets_graph[gdf_streets_graph[self._street_graph_fid] == adj_fid].iloc[0]

                # distance = gdf_streets_graph[gdf_streets_graph[self._street_graph_fid] == street_graph_fid].iloc[0]['geometry_prj'].distance(street_graph_adj['geometry_prj'])   # DOC: pure distance between two street geometries
                distance = street_distances[s_fid] + street_graph['geometry_prj'].centroid.distance(street_graph_adj['geometry_prj'].centroid) # DOC: distance between centroids of streets in path
                
                if distance <= max_distance:
                    streets_graph_fids.append(adj_fid)
                    street_distances[adj_fid] = distance
                
            del street_distances[s_fid]   # DOC: avoid considering distance to the same street again

        street_gates = street_gates_df.sort_values(by='distance').drop_duplicates(subset=self._gate_fid, keep='first').reset_index(drop=True).to_dict(orient='records') if not street_gates_df.empty else []
        return street_gates


    def __call__(self, **kwargs):

        # DOC: Extract kwargs values ----------------------------------------------------------------------------------
        gdf_wd = kwargs['gdf_water_depth']
        bbox = kwargs['bbox']
        t_srs = kwargs['t_srs']

        # DOC: Get street graph and street layer with which compute flood ---------------------------------------------
        gdf_streets_graph = module_retriever.retrieve_venezia_wfs(
            provider = self._street_graph_layer_id,
            bbox = bbox,
            buffer_points = False
        )        
        gdf_streets = module_retriever.retrieve_venezia_wfs(
            provider = self.street_layer_id,
            bbox = bbox,
            buffer_points = False
        )
        gdf_gate = module_retriever.retrieve_venezia_wfs(
            provider = self._gates_layer_id,
            bbox = bbox,
            buffer_points = False
        )
        gdf_gate.to_file('gdf_gate.gpkg', driver='GPKG')

        # DOC: Build G(V,E) from street graph layer --------------------------------------------------------------------
        """
        L'idea: le geometrie non sono un grafo e l'intersezione delle geometrie non funziona per la sua costruzione causa ponti e gallerie.
        Noto però che le strade hanno geometrie diverse per ogni incrocio:
        quindi per ogni geometria salvo punto iniziale e finale (endpoint), esplodo (2 record a geometria), self-join su endpoint.
        """
        gdf_nodes = gpd.GeoDataFrame(geometry = gdf_streets_graph.geometry.apply(lambda g: MultiPoint(points=[ list(g.geoms)[0].coords[0], list(g.geoms)[0].coords[-1] ])), crs=gdf_streets_graph.crs)
        gdf_nodes[self._street_graph_fid] = gdf_streets_graph[self._street_graph_fid]
        point_to_id = lambda p: int(''.join(map(str, list(*p.coords))).replace('.', ''))
        gdf_nodes = gdf_nodes.explode(ignore_index=True).reset_index(drop=True)
        gdf_nodes['endpoint'] = gdf_nodes.geometry.apply(lambda g: point_to_id(g))

        gdf_nodes = gdf_nodes.merge(gdf_nodes, left_on='endpoint', right_on='endpoint', suffixes=('_start', '_end'))
        gdf_nodes = gdf_nodes[gdf_nodes[f'{self._street_graph_fid}_start'] != gdf_nodes[f'{self._street_graph_fid}_end']].reset_index(drop=True)

        gdf_nodes[f'endpoint_{self._street_graph_fid}'] = gdf_nodes.apply(lambda r: [r[f'{self._street_graph_fid}_start'], r[f'{self._street_graph_fid}_end']], axis=1)
        gdf_nodes = gdf_nodes.explode(f'endpoint_{self._street_graph_fid}', ignore_index=True).reset_index(drop=True)
        gdf_nodes = gdf_nodes.groupby(f'endpoint_{self._street_graph_fid}').aggregate({f'{self._street_graph_fid}_start': list, f'{self._street_graph_fid}_end': list}).reset_index()
        gdf_nodes[f'adjacent_{self._street_graph_fid}'] = gdf_nodes.apply(lambda r: list(set(r[f'{self._street_graph_fid}_start'] + r[f'{self._street_graph_fid}_end'])), axis=1)
        gdf_nodes = gdf_nodes.drop(columns=[f'{self._street_graph_fid}_start', f'{self._street_graph_fid}_end']).drop_duplicates(subset=[f'endpoint_{self._street_graph_fid}'], ignore_index=True)
        gdf_streets_graph = gdf_streets_graph.merge(gdf_nodes, left_on=self._street_graph_fid, right_on=f'endpoint_{self._street_graph_fid}', suffixes=('', '_node'), how='left')
        gdf_streets_graph[f'adjacent_{self._street_graph_fid}'] = gdf_streets_graph[f'adjacent_{self._street_graph_fid}'].apply(lambda l: l if type(l) is list else [])  # DOC: ensure adjacent is a list
        
        # DOC: Create UTMxx projected CRS for spatial operations ------------------------------------------------------
        gdf_wd['geometry_prj'] = gdf_wd['geometry'].to_crs(_consts._EPSG_UTMxx)
        gdf_gate['geometry_prj'] = gdf_gate['geometry'].to_crs(_consts._EPSG_UTMxx)
        gdf_streets_graph['geometry_prj'] = gdf_streets_graph['geometry'].to_crs(_consts._EPSG_UTMxx).buffer(self.street_buffer)
        gdf_streets['geometry_prj'] = gdf_streets['geometry'].to_crs(_consts._EPSG_UTMxx).buffer(self.street_buffer)

        # DOC: Compute flooded streets --------------------------------------------------------------------------------
        wd_tree = STRtree(gdf_wd['geometry_prj'].values)
        street_in_wd = pd.DataFrame(wd_tree.query(gdf_streets['geometry_prj'].values, predicate='intersects').T, columns=['street_idx', 'wd_idx'])
        street_in_wd['flood_area'] = street_in_wd.apply(lambda r: gdf_streets.loc[r['street_idx']]['geometry_prj'].intersection(gdf_wd.loc[r['wd_idx']]['geometry_prj']).area, axis=1)
        street_in_wd = street_in_wd[street_in_wd['flood_area'] >= self.flood_area_thresh]
        gdf_streets['is_flooded'] = gdf_streets.index.isin(street_in_wd['street_idx'])
        Logger.debug(f"Number of flooded streets: {gdf_streets['is_flooded'].sum()} out of {len(gdf_streets)}")

        # DOC: Search gates for flooded streets -----------------------------------------------------------------------
        
        # DOC: Assign gate_idx to streets in the street graph
        gate_tree = STRtree(gdf_gate['geometry_prj'].buffer(self.gates_buffer).values)
        streets_graph_with_gate = pd.DataFrame(gate_tree.query(gdf_streets_graph['geometry_prj'].values, predicate='intersects').T, columns=['street_idx', 'gate_idx']).groupby('street_idx').agg(list).reset_index()
        gdf_streets_graph['gate_idx'] = [list() for _ in range(len(gdf_streets_graph))]
        gdf_streets_graph.loc[streets_graph_with_gate['street_idx'], 'gate_idx'] = pd.Series(streets_graph_with_gate['gate_idx'].values, index=streets_graph_with_gate['street_idx']).values

        # DOC: Foreach flooded street, search for gates in the street graph
        gdf_streets['gates_guard'] = [list() for _ in range(len(gdf_streets))]
        gdf_streets.loc[gdf_streets.is_flooded, 'gates_guard'] = gdf_streets[gdf_streets['is_flooded']].apply(
            lambda r: self.bfs_streets_gates(
                gdf_streets_graph = gdf_streets_graph,
                street_graph_fid = r[self._street_graph_fid],
                gdf_gate = gdf_gate,
                max_distance = self.max_distance
            ),
            axis=1
        )
        Logger.debug(f"## Found {gdf_streets['gates_guard'].apply(lambda g: len(g) > 0).sum()} flooded streets with gates out of {gdf_streets.is_flooded.sum()} total flooded streets.")

        # DOC: Prepare features collection to be returned -------------------------------------------------------------
        all_gates_found = gdf_streets['gates_guard'].explode().dropna().reset_index(drop=True).apply(lambda g: g[self._gate_fid]).unique().tolist()
        gdf_gates_guard = gdf_gate[gdf_gate[self._gate_fid].isin(all_gates_found)].reset_index(drop=True)
        gdf_gates_guard = _utils.safe_json_df(gdf_gates_guard.drop(columns=['geometry_prj']).to_crs(crs=t_srs))
        gdf_streets = _utils.safe_json_df(gdf_streets.drop(columns=['geometry_prj']).to_crs(crs=t_srs))
        gates_guard_collection = _utils.set_crs_feature_collection(gdf_gates_guard.to_geo_dict(), t_srs)
        streets_collection = _utils.set_crs_feature_collection(gdf_streets.to_geo_dict(), t_srs)

        return streets_collection, gates_guard_collection



_ADD_OPS = {

    NearbyPumps.name: {
        'class': NearbyPumps,
        'args': NearbyPumps.args,
        'providers': [
            _consts._VENEZIA_WFS_PROVIDER,
            _consts._VENEZIA_WFS_CRITICAL_SITES_PROVIDER,
            _consts._OVERTURE_PROVIDER
        ]
    },

    AlertMethod.name: {
        'class': AlertMethod,
        'args': AlertMethod.args,
        'providers': [
            _consts._VENEZIA_WFS_PROVIDER,
            _consts._VENEZIA_WFS_CRITICAL_SITES_PROVIDER,
            _consts._OVERTURE_PROVIDER
        ]
    },

    GatesGuard.name: {
        'class': GatesGuard,
        'args': GatesGuard.args,
        'providers': [
            _consts._VENEZIA_WFS_PROVIDER,
            _consts._VENEZIA_WFS_CRITICAL_SITES_PROVIDER,
            _consts._OVERTURE_PROVIDER
        ]
    }
}