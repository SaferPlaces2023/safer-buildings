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

    name = 'nearby_pumps'
    args = [
         # DOC: The buffer to use around water depth areas (in meters). Default is 100.0 meters.
        'wd_buffer',
        'max_distance'          
    ]

    _pumps_layer_id = 'mv_risorse_p0109103_pompe'
    _pumps_identifier = 'gid'
    _nearby_pumps_basic_attributes = ['id', '_fid', 'gid', 'modello', 'indirizzo', 'longitude', 'latitude']

    def __init__(self, wd_buffer: float = 100.0, max_distance: float = 100.0):
        super().__init__(name=self.name)
        self._configure(wd_buffer=wd_buffer, max_distance=max_distance)

    def _configure(self, wd_buffer: float = 100.0, max_distance: float = 100.0):
        self.wd_buffer = float(wd_buffer) if type(wd_buffer) is str else wd_buffer
        self.max_distance = float(max_distance) if type(max_distance) is str else max_distance


    def opt(self, **kwargs):

        # DOC: Extract kwargs values
        gdf_buildings = kwargs['gdf_buildings']
        gdf_wd = kwargs['gdf_water_depth']
        t_srs = kwargs['t_srs']

        # DOC: Retrieve the pumps from the WFS provider
        gdf_pumps = module_retriever.retrieve_venezia_wfs(
            provider = f'{_consts._VENEZIA_WFS_PROVIDER}/{self._pumps_layer_id}',
            bbox = kwargs['bbox'],
            buffer_points = False 
        ).reset_index(drop=True)
        # DOC: Prepare a projected-crs geometry column to use in spatial operations
        gdf_pumps['geometry_3857'] = gdf_pumps.to_crs(epsg=3857).geometry
        gdf_pumps['geometry'] = gdf_pumps.geometry_3857.centroid
        # DOC: Convert the pumps GeoDataFrame to the target CRS
        gdf_pumps = gdf_pumps.to_crs(t_srs)
        gdf_pumps['longitude'] = gdf_pumps.geometry.x
        gdf_pumps['latitude'] = gdf_pumps.geometry.y

        # DOC: Convert CRS to Projected CRS (EPSG:3857) due to calculating distances
        gdf_buildings['geometry_3857'] = gdf_buildings.to_crs(epsg=3857).geometry

        # DOC: Create a spatial index for the water depth polygons
        wd_3857_tree = STRtree(gdf_wd.to_crs(epsg=3857).geometry.values)

        # DOC: Find pumps whithin a distance from the water depth polygons
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
        Logger.debug(f"## Found {len(gdf_pumps)} nearby pumps for the flooded areas.")

        # DOC: Find wd area for each building
        buildings_in_wd_query = wd_3857_tree.query(gdf_buildings['geometry_3857'], predicate='intersects', distance=0)
        bld_wd_idxs = pd.DataFrame(buildings_in_wd_query.T, columns=['building_idx', 'wd_idxs']).groupby('building_idx').agg(list).reset_index()

        # region: con join
        gdf_buildings = pd.merge(
            gdf_buildings, bld_wd_idxs,
            left_index=True, right_on='building_idx', how='left'
        ).drop(columns='building_idx')
        gdf_buildings_pumps = gdf_buildings.explode('wd_idxs')
        gdf_buildings_pumps = pd.merge(
            gdf_buildings_pumps, gdf_pumps[self._nearby_pumps_basic_attributes + ['geometry_3857']],
            left_on='wd_idxs', right_index=True, how='left',
            suffixes=('', '_pump')
        ).drop(columns='wd_idxs')
        gdf_buildings_pumps = gdf_buildings_pumps[gdf_buildings_pumps['geometry_3857'].notna()]
        Logger.debug('##### joined')
        gdf_buildings_pumps['distance_pump'] = distance(gdf_buildings_pumps['geometry_3857'].values, gdf_buildings_pumps['geometry_3857_pump'].values)
        Logger.debug('##### distance calculated')
        gdf_buildings_pumps = gdf_buildings_pumps[gdf_buildings_pumps['distance_pump'] <= self.max_distance]
        gdf_buildings_pumps = gdf_buildings_pumps.groupby(gdf_buildings_pumps.index).agg({
            col: 'first' for col in gdf_buildings_pumps.columns if not col.endswith('_pump')
        } | {
            col: list for col in gdf_buildings_pumps.columns if col.endswith('_pump')
        }).rename(columns={
            col: col.replace('_pump', '') for col in gdf_buildings_pumps.columns if col.endswith('_pump')
        }).drop(columns=['geometry_3857'])
        buildings_nearby_pumps = gdf_buildings_pumps.to_dict(orient='records')
        gdf_buildings['nearby_pumps'] = [list() for _ in range(len(gdf_buildings))]
        gdf_buildings.loc[gdf_buildings_pumps.index, 'nearby_pumps'] = pd.Series(buildings_nearby_pumps, index=gdf_buildings_pumps.index)


        # endregion: con join

        # bld_wd_idxs.set_index('building_idx', inplace=True)
        # # DOC: For each building, assign the nearby pumps selected for the building water depth area
        # def get_building_pumps(building):
        #     if building.Index not in bld_wd_idxs.index:
        #         return []
        #     wd_idxs = list(set(bld_wd_idxs.loc[building.Index]['wd_idxs']) & set(gdf_pumps.index))
        #     if len(wd_idxs) < 1:
        #         return []
        #     bld_pumps = gdf_pumps.loc[wd_idxs]
        #     if bld_pumps.empty:
        #         return []
        #     bld_pumps['distance'] = bld_pumps['geometry_3857'].distance(building.geometry_3857)
        #     bld_pumps = bld_pumps[bld_pumps['distance'] <= self.max_distance]
        #     if bld_pumps.empty:
        #         return []
        #     bld_pumps = bld_pumps[self._nearby_pumps_basic_attributes].to_dict(orient='records')
        #     return bld_pumps
        # gdf_buildings['nearby_pumps'] = [
        #     get_building_pumps(buidling)
        #     for buidling in tqdm(
        #         gdf_buildings.itertuples(), 
        #         total=len(gdf_buildings), mininterval=0.5, disable=not is_debug_mode(), desc="Assign nearby pumps to each building"
        #     )
        # ]

        # DOC: Convert the nearby pumps GeoDataFrame to the target CRS
        gdf_pumps = gdf_pumps.reset_index(drop=True).drop(columns=['geometry_3857'])
        gdf_pumps = _utils.df_dt_col_to_isoformat(gdf_pumps)
        nearby_pumps_collection = gdf_pumps.to_geo_dict()
        Logger.debug(f"## Buildings with nearby pumps: {gdf_buildings[gdf_buildings['is_flooded']]['nearby_pumps'].apply(lambda npumps: len(npumps)>0).sum()} out of {len(gdf_buildings[gdf_buildings['is_flooded']])} total flooded buildings.")

        # DOC: Return the updated GeoDataFrame with the nearest pumps and the nearby pumps collection
        return gdf_buildings, nearby_pumps_collection





    def __call__(self, **kwargs):

        return self.opt(**kwargs)

        gdf_buildings = kwargs['gdf_buildings']
        gdf_wd = kwargs['gdf_water_depth']

        # DOC: Retrieve the pumps from the WFS provider
        gdf_pumps = module_retriever.retrieve_venezia_wfs(
            provider = f'{_consts._VENEZIA_WFS_PROVIDER}/{self._pumps_layer_id}',
            bbox = kwargs['bbox'],
            buffer_points=False 
        )

        # DOC: Convert CRS to Projected CRS (EPSG:3857) due to calculating distances
        gdf_buildings_3857 = gdf_buildings.to_crs(epsg=3857)
        gdf_wd_3857 = gdf_wd.to_crs(epsg=3857)
        if self.wd_buffer > 0:
            gdf_wd_3857['geometry'] = gdf_wd_3857.geometry.buffer(self.wd_buffer)
        gdf_pumps_3857 = gdf_pumps.to_crs(epsg=3857)
        gdf_pumps_3857['geometry'] = gdf_pumps_3857.centroid
        # DOC: Convert the pumps GeoDataFrame to the target CRS
        gdf_pumps_tsrs  = gdf_pumps.to_crs(kwargs['t_srs'])

        # DOC: Build a search tree for pumps geometries        
        pumps_3857_tree = STRtree(gdf_pumps_3857.geometry.values)
        pump_geom_to_index = {id(g): i for i, g in enumerate(gdf_pumps_3857.geometry.values)}
        def pumps_in_flood_area(wd_geom):
            candidates = pumps_3857_tree.geometries.take(pumps_3857_tree.query(wd_geom)).tolist()
            return gdf_pumps_3857.iloc[[pump_geom_to_index[id(g)] for g in candidates if wd_geom.intersects(g)]]

        # DOC: Foreach water depth area, retrieve (if any) the pump within the water depth area
        Logger.debug(f"### Retrieving nearby pumps for {len(gdf_wd_3857)} water depth areas.")
        nearby_pumps_identifiers = set()
        gdf_buildings['nearby_pumps'] = [list() for _ in range(len(gdf_buildings))]
        # for i_wd, wd in gdf_wd_3857.iterrows():
        for i_wd, wd in tqdm(gdf_wd_3857.iterrows(), total=len(gdf_wd_3857), mininterval=0.5, disable=not is_debug_mode(), desc="Search for each flooded area"):

            # DOC: Get the pumps within the water depth area
            wd_pumps = pumps_in_flood_area(wd.geometry)
            if wd_pumps.empty:
                continue
            wd_pumps = wd_pumps.drop_duplicates(subset=self._pumps_identifier)

            # DOC: Collect the nearby pumps GeoDataFrame in the target CRS
            wd_pumps = gdf_pumps_tsrs.loc[wd_pumps.index]
            nearby_pumps_identifiers.update(wd_pumps[self._pumps_identifier].to_list())

            # DOC: Mark buildings in the water depth area affected by retrieved nearby pumps
            affected_buildings = gdf_buildings_3857[gdf_buildings_3857.geometry.intersects(wd.geometry)]
            if affected_buildings.empty:
                continue
            wd_pumps['location'] = wd_pumps['geometry'].apply(lambda geom: [geom.centroid.x, geom.centroid.y])
            wd_pumps = pd.DataFrame(wd_pumps.drop(columns='geometry'))
            wd_pumps = wd_pumps[self._nearby_pumps_basic_attributes + ['location']]
            wd_pumps = _utils.df_dt_col_to_isoformat(wd_pumps)
            wd_pumps = json.loads(wd_pumps.to_json(orient='records'))
            affected_building_index = affected_buildings.index.to_list()
            gdf_buildings.loc[affected_building_index, 'nearby_pumps'] = pd.Series([wd_pumps for _ in range(len(affected_building_index))], index=affected_building_index)

        # DOC: If nearby pumps were found, build a feature collection of them
        if len(nearby_pumps_identifiers) > 0:
            nearby_pumps_gdf = gdf_pumps_tsrs[gdf_pumps_tsrs[self._pumps_identifier].isin(nearby_pumps_identifiers)]
            nearby_pumps_gdf = _utils.df_dt_col_to_isoformat(nearby_pumps_gdf)
            nearby_pumps_collection = nearby_pumps_gdf.to_geo_dict()
        else:
            nearby_pumps_collection = { 'type': 'FeatureCollection', 'features': [] }

        Logger.debug(f"## Found {len(nearby_pumps_collection['features'])} nearby pumps for the flooded areas.")
        
        # DOC: Return the updated GeoDataFrame with the nearest pumps and the nearby pumps collection
        return gdf_buildings, nearby_pumps_collection

            



# DOC: This is for VENEZIA_WFS_PROVIDER (Alert method to use)
class AlertMethod(AdditionalOperation):

    name = 'alert_method'
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
        
        gdf_buildings = kwargs['gdf_buildings']
        gdf_wd = kwargs['gdf_water_depth']

        # DOC: Retrieve the alert area and method (only related to flooding event) from the WFS provider
        gdf_alert_area = module_retriever.retrieve_venezia_wfs(
            provider = f'{_consts._VENEZIA_WFS_PROVIDER}/{self._alert_area_layer_id}',
            bbox = kwargs['bbox']
        )
        gdf_alert_method = module_retriever.retrieve_venezia_wfs(
            provider = f'{_consts._VENEZIA_WFS_PROVIDER}/{self._alert_method_layer_id}',
            bbox = kwargs['bbox']
        )
        gdf_alert_method = gdf_alert_method[(gdf_alert_method['allagament_t']=='Sì') | (gdf_alert_method['r_altro']=='Acqua Alta')]
        # DOC: Convert the alert mtehod GeoDataFrame to the target CRS
        gdf_aert_method_tsrs = gdf_alert_method.to_crs(kwargs['t_srs'])

        # DOC: Convert CRS to Projected CRS (EPSG:3857) due to geometry calculations
        gdf_buildings_3857 = gdf_buildings.to_crs(epsg=3857)
        gdf_wd_3857 = gdf_wd.to_crs(epsg=3857)
        if self.wd_buffer > 0:
            gdf_wd_3857['geometry'] = gdf_wd_3857.geometry.buffer(self.wd_buffer)
        gdf_alert_area_3857 = gdf_alert_area.to_crs(epsg=3857)
        gdf_alert_method_3857 = gdf_alert_method.to_crs(epsg=3857)
        gdf_alert_method_3857['geometry'] = gdf_alert_method_3857.centroid

        # DOC: Build a search tree for buildings geometries
        buildings_3857_tree = STRtree(gdf_buildings_3857.geometry.values)
        building_geom_to_index = {id(g): i for i, g in enumerate(gdf_buildings_3857.geometry.values)}
        def buildings_in_alert_area(alert_area_geom):
            candidates = buildings_3857_tree.geometries.take(buildings_3857_tree.query(alert_area_geom)).tolist()
            return gdf_buildings_3857.iloc[[building_geom_to_index[id(g)] for g in candidates if alert_area_geom.intersects(g)]]

        # DOC: Build a search tree for alert area geometries
        alert_area_3857_tree = STRtree(gdf_alert_area_3857.geometry.values)
        alert_area_geom_to_index = {id(g): i for i, g in enumerate(gdf_alert_area_3857.geometry.values)}
        def alert_area_in_flood_area(wd_geom):
            candidates = alert_area_3857_tree.geometries.take(alert_area_3857_tree.query(wd_geom)).tolist()
            return gdf_alert_area_3857.iloc[[alert_area_geom_to_index[id(g)] for g in candidates if wd_geom.intersects(g)]]
        
        # DOC: Build a search tree for alert method geometries
        alert_method_3857_tree = STRtree(gdf_alert_method_3857.geometry.values)
        alert_method_geom_to_index = {id(g): i for i, g in enumerate(gdf_alert_method_3857.geometry.values)}
        def alert_method_in_alert_area(alert_area_geom):
            candidates = alert_method_3857_tree.geometries.take(alert_method_3857_tree.query(alert_area_geom)).tolist()
            return gdf_alert_method_3857.iloc[[alert_method_geom_to_index[id(g)] for g in candidates if alert_area_geom.intersects(g)]]

        # DOC: Foreach flood area, get the relative alert method
        Logger.debug(f"### Retrieving alert methods for {len(gdf_wd_3857)} water depth areas.")
        alert_method_identifiers = set()
        gdf_buildings['alert_method'] = [list() for _ in range(len(gdf_buildings))]
        # for i_wd, wd in gdf_wd_3857.iterrows():
        for i_wd, wd in tqdm(gdf_wd_3857.iterrows(), total=len(gdf_wd_3857), mininterval=0.5, disable=not is_debug_mode(), desc="Search for each flooded area"):
            
            # DOC: Get the alert area within the water depth area
            wd_alert_area = alert_area_in_flood_area(wd.geometry)
            if wd_alert_area.empty:
                continue
            wd_alert_area = wd_alert_area.drop_duplicates(subset=self._alert_area_identifier)

            # DOC: Get the alert method within the water depth area
            wd_alert_method = alert_method_in_alert_area(wd.geometry)   # !!!: There is no an explicit relation between alert area and method, so we look at which method is contained in the alert area
            if wd_alert_method.empty:
                continue    # DOC: this should only when area is related to a not-flooding alert method or if method is not acoustic

            wd_alert_method = wd_alert_method.drop_duplicates(subset=self._alert_method_identifier)
            
            # DOC: Collect the alert method GeoDataFrame
            wd_alert_method = gdf_aert_method_tsrs.loc[wd_alert_method.index]
            alert_method_identifiers.update(wd_alert_method[self._alert_method_identifier].to_list())

            # DOC: Mark (flooded) buildings in the alert area with relative alert method
            alert_buildings = buildings_in_alert_area(wd_alert_area.geometry.union_all())
            # ???: Filter only alert_buildings.is_flooded → (logically yes but maybe it is important to have the alert method for all buildings in the alert area)
            if alert_buildings.empty:
                continue
            wd_alert_method['location'] = wd_alert_method['geometry'].apply(lambda geom: [geom.centroid.x, geom.centroid.y])
            wd_alert_method = pd.DataFrame(wd_alert_method.drop(columns='geometry'))
            wd_alert_method = wd_alert_method[self._alert_method_basic_attributes + ['location']]
            wd_alert_method = _utils.df_dt_col_to_isoformat(wd_alert_method)
            wd_alert_method = json.loads(wd_alert_method.to_json(orient='records'))
            alert_building_index = alert_buildings.index.to_list()
            gdf_buildings.loc[alert_building_index, 'alert_method'] = pd.Series([wd_alert_method for _ in range(len(alert_building_index))], index=alert_building_index)

        # DOC: If alert methods weere founde, build a feature colection of them
        if len(alert_method_identifiers) > 0:
            alert_method_gdf = gdf_aert_method_tsrs[gdf_aert_method_tsrs[self._alert_method_identifier].isin(alert_method_identifiers)]
            alert_method_gdf = _utils.df_dt_col_to_isoformat(alert_method_gdf)
            alert_method_collection = alert_method_gdf.to_geo_dict()
        else:
            alert_method_collection = { 'type': 'FeatureCollection', 'features': [] }

        Logger.debug(f"## Found {len(alert_method_collection['features'])} alert methods for the flooded areas.")

        # DOC: Return the updated building GeoDataFrame with the alert method and the alert method collection    
        return gdf_buildings, alert_method_collection





_ADD_OPS = {

    NearbyPumps.name: {
        'class': NearbyPumps,
        'args': NearbyPumps.args,
        'providers': [
            _consts._VENEZIA_WFS_PROVIDER,
            _consts._OVERTURE_PROVIDER
        ]
    },

    AlertMethod.name: {
        'class': AlertMethod,
        'args': AlertMethod.args,
        'providers': [
            _consts._VENEZIA_WFS_PROVIDER,
            _consts._OVERTURE_PROVIDER
        ]
    }
}