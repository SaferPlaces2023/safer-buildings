import os
import json
import datetime
from tqdm import tqdm

import pandas as pd
import geopandas as gpd
from shapely import distance
from shapely.geometry import box, MultiPolygon
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
        
        # DOC: Convert CRS to Projected CRS (EPSG:3857) due to calculating distances
        gdf_buildings['geometry_3857'] = gdf_buildings.to_crs(epsg=3857).geometry

        # DOC: Retrieve the pumps from the WFS provider ---------------------------------------------------------------
        gdf_pumps = module_retriever.retrieve_venezia_wfs(
            provider = f'{_consts._VENEZIA_WFS_PROVIDER}/{self._pumps_layer_id}',
            bbox = bbox,
            buffer_points = False 
        ).reset_index(drop=True)
        # DOC: Prepare a projected-crs geometry column to use in spatial operations
        gdf_pumps['geometry_3857'] = gdf_pumps.to_crs(epsg=3857).geometry
        gdf_pumps['geometry'] = gdf_pumps['geometry_3857'].centroid
        # DOC: Convert the pumps GeoDataFrame to the target CRS
        gdf_pumps = gdf_pumps.to_crs(t_srs)
        gdf_pumps['longitude'] = gdf_pumps.geometry.x
        gdf_pumps['latitude'] = gdf_pumps.geometry.y
        self._nearby_pumps_basic_attributes.extend(['longitude', 'latitude'])

        # DOC: Create a spatial index for the water depth polygons ----------------------------------------------------
        wd_3857_tree = STRtree(gdf_wd.to_crs(epsg=3857).geometry.values)

        # DOC: Find pumps whithin a distance from the water depth polygons --------------------------------------------
        pumps_in_wd_query = wd_3857_tree.query(gdf_pumps['geometry_3857'], predicate='dwithin', distance=self.max_distance)
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
        buildings_in_wd_query = wd_3857_tree.query(gdf_buildings['geometry_3857'], predicate='intersects', distance=0)
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
            gdf_buildings_pumps, gdf_pumps[[f'{attr}_pump' for attr in self._nearby_pumps_basic_attributes] + ['geometry_3857_pump']],
            left_on='wd_idxs', right_index=True, how='left'
        ).drop(columns='wd_idxs')
        gdf_buildings_pumps = gdf_buildings_pumps[gdf_buildings_pumps['geometry_3857_pump'].notna()]
        
        # DOC: Calculate the distance between each building and its related pumps (then filter by max_distance)
        gdf_buildings_pumps['distance_pump'] = distance(gdf_buildings_pumps['geometry_3857'].values, gdf_buildings_pumps['geometry_3857_pump'].values)
        gdf_buildings_pumps = gdf_buildings_pumps[gdf_buildings_pumps['distance_pump'] <= self.max_distance]
        
        # DOC: Foreach building, collect the nearby pumps (only basic attributes) in a list
        gdf_buildings_pumps = gdf_buildings_pumps[[col for col in gdf_buildings_pumps.columns if col.endswith('_pump')]].drop(columns=['geometry_3857_pump'])
        gdf_buildings_pumps.rename(columns={col: col.replace('_pump', '') for col in gdf_buildings_pumps.columns}, inplace=True)
        gdf_buildings_pumps = gdf_buildings_pumps.groupby(gdf_buildings_pumps.index).agg(list)
        buildings_nearby_pumps = gdf_buildings_pumps.apply(lambda row: pd.DataFrame(row.to_dict()).to_dict(orient='records'), axis=1).to_list()
        # DOC: Add the nearby pumps to the buildings GeoDataFrame       
        gdf_buildings['nearby_pumps'] = [list() for _ in range(len(gdf_buildings))]
        gdf_buildings.loc[gdf_buildings_pumps.index, 'nearby_pumps'] = pd.Series(buildings_nearby_pumps, index=gdf_buildings_pumps.index)
        gdf_buildings.drop(columns=['geometry_3857'], inplace=True)

        # DOC: Prepare the pumps GeoDataFrame to be returned ----------------------------------------------------------
        gdf_pumps.rename(columns={col: col.replace('_pump', '') for col in gdf_pumps.columns}, inplace=True)
        gdf_pumps.drop_duplicates(subset=self._pumps_identifier, inplace=True)
        gdf_pumps = gdf_pumps.drop(columns=['geometry_3857']).reset_index(drop=True)
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

        # DOC: Convert CRS to Projected CRS (EPSG:3857) due to compute spatial operations
        gdf_buildings['geometry_3857'] = gdf_buildings.to_crs(epsg=3857).geometry

        # DOC: Retrieve the alert areas -------------------------------------------------------------------------------
        gdf_alert_area = module_retriever.retrieve_venezia_wfs(
            provider = f'{_consts._VENEZIA_WFS_PROVIDER}/{self._alert_area_layer_id}',
            bbox = bbox,
            buffer_points=False 
        ).reset_index(drop=True)
        gdf_alert_area['geometry_3857'] = gdf_alert_area.to_crs(epsg=3857).geometry
        # DOC: Retrieve the alert methods (only related to flooding event) from the WFS provider
        gdf_alert_method = module_retriever.retrieve_venezia_wfs(
            provider = f'{_consts._VENEZIA_WFS_PROVIDER}/{self._alert_method_layer_id}',
            bbox = bbox,
            buffer_points = False
        )
        gdf_alert_method = gdf_alert_method[(gdf_alert_method['allagament_t']=='Sì') | (gdf_alert_method['r_altro']=='Acqua Alta')].reset_index(drop=True)
        gdf_alert_method['geometry_3857'] = gdf_alert_method.to_crs(epsg=3857).geometry
        # DOC: Convert the alert method GeoDataFrame to the target CRS
        gdf_alert_method['geometry'] = gpd.GeoSeries(gdf_alert_method['geometry_3857'].centroid, crs='EPSG:3857').to_crs(gdf_alert_method.crs)
        gdf_alert_method = gdf_alert_method.to_crs(t_srs)
        gdf_alert_method['longitude'] = gdf_alert_method.geometry.x
        gdf_alert_method['latitude'] = gdf_alert_method.geometry.y
        self._alert_method_basic_attributes.extend(['longitude', 'latitude'])

        # DOC: Filter out alert area that not contains alert methods
        gdf_alert_area = gdf_alert_area[gdf_alert_area['geometry_3857'].intersects(gdf_alert_method['geometry_3857'].union_all())].reset_index(drop=True)
        
        # DOC: Create a spatial index for the water depth polygons ----------------------------------------------------
        wd_tree = STRtree(gdf_wd.to_crs(epsg=3857).buffer(max(self.wd_buffer, 0)).geometry.values)

        # DOC: Alert area in water depth areas ------------------------------------------------------------------------
        alert_area_in_wd_query = wd_tree.query(gdf_alert_area['geometry_3857'].values, predicate='intersects', distance=0)
        alert_area_selection = gdf_alert_area.iloc[alert_area_in_wd_query[0,:]]
        alert_area_selection['index'] = alert_area_selection.index
        alert_area_selection.drop_duplicates(subset='index', inplace=True)
        alert_area_selection.drop(columns='index', inplace=True)
        Logger.debug(f"## Found {len(alert_area_selection)} alert areas in the water depth areas.")

        # DOC: Alert method in alert area -----------------------------------------------------------------------------
        alert_area_selection_tree = STRtree(alert_area_selection['geometry_3857'].values)
        alert_method_query = alert_area_selection_tree.query(gdf_alert_method['geometry_3857'].values, predicate='within')
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
        buildings_alert_query = alert_area_selection_tree.query(gdf_buildings['geometry_3857'].values, predicate='intersects')
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
            gdf_buildings_alert, gdf_alert_method[[f'{attr}_alert_method' for attr in self._alert_method_basic_attributes] + ['geometry_3857_alert_method']],
            left_on='alert_area_idx', right_index=True, how='left'
        ).drop(columns='alert_area_idx')
        gdf_buildings_alert = gdf_buildings_alert[gdf_buildings_alert['geometry_3857_alert_method'].notna()]

        # DOC: Foreach building, collect the alert methods (only basic attributes) in a list
        gdf_buildings_alert = gdf_buildings_alert[[col for col in gdf_buildings_alert.columns if col.endswith('_alert_method')]].drop(columns=['geometry_3857_alert_method'])
        gdf_buildings_alert.rename(columns={col: col.replace('_alert_method', '') for col in gdf_buildings_alert.columns}, inplace=True)
        gdf_buildings_alert = gdf_buildings_alert.groupby(gdf_buildings_alert.index).agg(list)
        buildings_alert_methods = gdf_buildings_alert.apply(lambda row: pd.DataFrame(row.to_dict()).to_dict(orient='records'), axis=1).to_list()
        # DOC: Add the alert methods to the buildings GeoDataFrame
        gdf_buildings['alert_method'] = [list() for _ in range(len(gdf_buildings))]
        gdf_buildings.loc[gdf_buildings_alert.index, 'alert_method'] = pd.Series(buildings_alert_methods, index=gdf_buildings_alert.index)
        gdf_buildings.drop(columns=['geometry_3857'], inplace=True)

        # DOC: Prepare the alert method GeoDataFrame to be returned -------------------------------------------------
        gdf_alert_method.rename(columns={col: col.replace('_alert_method', '') for col in gdf_alert_method.columns}, inplace=True)
        gdf_alert_method.drop_duplicates(subset=self._alert_method_identifier, inplace=True)
        gdf_alert_method = gdf_alert_method.drop(columns=['geometry_3857', 'alert_method_idx']).reset_index(drop=True)
        gdf_alert_method = _utils.safe_json_df(gdf_alert_method)
        alert_method_collection = gdf_alert_method.to_geo_dict()
        alert_method_collection = _utils.set_crs_feature_collection(alert_method_collection, t_srs)
        Logger.debug(f"## Buildings with alert methods: {gdf_buildings[gdf_buildings[_consts._COL_IS_FLOODED]]['alert_method'].apply(lambda am: len(am)>0).sum()} out of {len(gdf_buildings[gdf_buildings[_consts._COL_IS_FLOODED]])} total flooded buildings.")

        # DOC: Return the updated GeoDataFrame with the alert method and the alert method collection ----------------
        return gdf_buildings, alert_method_collection





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
    }
}