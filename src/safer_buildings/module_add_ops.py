import os
import json
import datetime

import pandas as pd
import geopandas as gpd
from shapely.geometry import box

from . import _consts, _utils
from . import module_retriever
from .module_log import Logger



class AdditionalOperation():

    def __init__(self, name: str):
        self.name = name

    def __str__(self):
        return f"Callable: {self.name}"
        


# DOC: This is for VENEZIA_WFS_PROVIDER (Find the nearest pump to a flooded building)
class NearbyPumps(AdditionalOperation):

    name = 'nearby-pumps'
    args = [
        # DOC: Maximum distance to consider a pump as "nearby" (in meters). Default is 1000.0 meters.
        'max_distance',
        # DOC: If True, only consider pumps that intersect with the water depth area of the building. Default is False.
        'wd_intersect'              
    ]

    _pumps_layer_id = 'mv_risorse_p0109103_pompe'
    _nearby_pumps_basic_attributes = ['id', 'modello', 'distance', 'indirizzo']

    def __init__(self, max_distance: float = 1000.0, wd_intersect: bool = False):
        super().__init__(name=self.name)
        self._configure(
            max_distance = max_distance, 
            wd_intersect = wd_intersect
        )
        

    def _configure(self, max_distance: float = 1000.0, wd_intersect: bool = False):
        self.max_distance = float(max_distance) if type(max_distance) is str else max_distance
        self.wd_intersect = bool(wd_intersect) if type(wd_intersect) is str else wd_intersect


    def __call__(self, **kwargs):

        gdf_buildings = kwargs['gdf_buildings']
        gdf_wd = kwargs['gdf_water_depth']

        # DOC: Retrieve the pumps from the WFS provider
        gdf_pumps = module_retriever.retrieve_venezia_wfs(
            provider = f'{_consts._VENEZIA_WFS_PROVIDER}/{self._pumps_layer_id}',
            bbox = kwargs['bbox']
        )

        # DOC: Convert CRS to Projected CRS (EPSG:3857) due to calculating distances
        gdf_buildings_3857 = gdf_buildings.to_crs(epsg=3857)
        gdf_wd_3857 = gdf_wd.to_crs(epsg=3857)
        gdf_pumps_3857 = gdf_pumps.to_crs(epsg=3857)
        gdf_pumps_3857['geometry'] = gdf_pumps_3857.centroid


        # DOC: Foreach building, retrieve (if any) the pump within the water depth area
        nearby_pumps_gdfs = []
        gdf_buildings['nearby_pumps'] = [list() for _ in range(len(gdf_buildings))]
        for i_building, building in gdf_buildings_3857.iterrows():

            # DOC: Skip if the building is not flooded
            if not building.is_flooded:
                continue
        
            if self.wd_intersect:
                # DOC: Get the flooded area for the building
                building_bbox = building.geometry.bounds
                building_wd = gdf_wd_3857.cx[building_bbox[0]:building_bbox[2], building_bbox[1]:building_bbox[3]]
                building_wd_bbox = building_wd.total_bounds
                # DOC: Get the pumps within the building's water depth area
                building_wd_pumps = gdf_pumps_3857.cx[building_wd_bbox[0]:building_wd_bbox[2], building_wd_bbox[1]:building_wd_bbox[3]]
            else:
                # DOC: Get the pumps within the building's buffer area bounding box
                building_bbox = building.geometry.buffer(self.max_distance).bounds
                building_wd_pumps = gdf_pumps_3857.cx[building_bbox[0]:building_bbox[2], building_bbox[1]:building_bbox[3]]
                
            if building_wd_pumps.empty:
                continue

            building_wd_pumps = building_wd_pumps.drop_duplicates(subset='id')

            # DOC: Calculate the distance from the building to each pump
            building_wd_pumps['distance'] = building.geometry.distance(building_wd_pumps.geometry)
            # DOC: Filter pumps within the max distance
            nearby_pumps = building_wd_pumps[building_wd_pumps['distance'] <= self.max_distance]
            if nearby_pumps.empty:
                continue
            
            # DOC: Collect the nearby pumps GeoDataFrame
            nearby_pumps = nearby_pumps.to_crs(kwargs['t_srs'])
            nearby_pumps_gdfs.append(nearby_pumps)
            
            # DOC: Add pump information to the building
            nearby_pumps['location'] = nearby_pumps['geometry'].apply(lambda geom: [geom.x, geom.y])
            nearby_pumps = pd.DataFrame(nearby_pumps.drop(columns='geometry'))
            nearby_pumps = nearby_pumps[self._nearby_pumps_basic_attributes + ['location']]
            nearby_pumps = _utils.df_dt_col_to_isoformat(nearby_pumps)
            nearby_pumps = json.loads(nearby_pumps.to_json(orient='records'))
            gdf_buildings.at[i_building, 'nearby_pumps'] = nearby_pumps


        if nearby_pumps_gdfs:
            nearby_pumps_gdf = pd.concat(nearby_pumps_gdfs, ignore_index=True).drop_duplicates(subset='id').reset_index(drop=True)
            nearby_pumps_gdf = _utils.df_dt_col_to_isoformat(nearby_pumps_gdf)
            nearby_pumps_collection = nearby_pumps_gdf.to_geo_dict()
        else:
            nearby_pumps_collection = { 'type': 'FeatureCollection', 'features': [] }

        # DOC: Return the updated GeoDataFrame with the nearest pumps and the nearby pumps collection
        return gdf_buildings, nearby_pumps_collection



# DOC: This is for VENEZIA_WFS_PROVIDER (Alert method to use)
class AlertMethod(AdditionalOperation):

    name = 'alert-method'
    args = [
        # DOC: The buffer to use around water depth areas (in meters). Default is 100.0 meters.
        'wd_buffer'
    ]

    _alert_area_layer_id = 'v_pc_p0103011_allertamento'
    _alert_method_layer_id = 'v_pc_p0103013_allertamento'
    _alert_method_basic_attributes = ['id', 'denom', 'indirizzo', 'strumento_t']

    def __init__(self, wd_buffer: float = 100.0):
        super().__init__(name=self.name)
        self._configure(wd_buffer=wd_buffer)

    def _configure(self, wd_buffer: float = 100.0):
        self.wd_buffer = float(wd_buffer) if type(wd_buffer) is str else wd_buffer

    def __call__(self, **kwargs):
        
        gdf_buildings = kwargs['gdf_buildings']
        gdf_wd = kwargs['gdf_water_depth']

        # DOC: Retieve the alert area and method (only related to flooding event) from the WFS provider
        gdf_alert_area = module_retriever.retrieve_venezia_wfs(
            provider = f'{_consts._VENEZIA_WFS_PROVIDER}/{self._alert_area_layer_id}',
            bbox = kwargs['bbox']
        )
        gdf_alert_method = module_retriever.retrieve_venezia_wfs(
            provider = f'{_consts._VENEZIA_WFS_PROVIDER}/{self._alert_method_layer_id}',
            bbox = kwargs['bbox']
        )
        gdf_alert_method = gdf_alert_method[(gdf_alert_method['allagament_t']=='Sì') | (gdf_alert_method['r_altro']=='Acqua Alta')]

        # DOC: Convert CRS to Projected CRS (EPSG:3857) due to calculating distances
        gdf_buildings_3857 = gdf_buildings.to_crs(epsg=3857)
        gdf_wd_3857 = gdf_wd.to_crs(epsg=3857)
        gdf_alert_area_3857 = gdf_alert_area.to_crs(epsg=3857)
        gdf_alert_method_3857 = gdf_alert_method.to_crs(epsg=3857)
        gdf_alert_method_3857['geometry'] = gdf_alert_method_3857.centroid

        # DOC: Foreach flood area, get the relative alert method
        alert_method_gdfs = []
        gdf_buildings['alert_method'] = [list() for _ in range(len(gdf_buildings))]
        for i_wd, wd in gdf_wd_3857.iterrows():
            
            # DOC: Get the alert area within the water depth area
            wd_bounds = wd.geometry.buffer(self.wd_buffer).bounds if self.wd_buffer > 0 else wd.geometry.bounds
            wd_alert_area = gdf_alert_area_3857.cx[wd_bounds[0]:wd_bounds[2], wd_bounds[1]:wd_bounds[3]]
            if wd_alert_area.empty:
                continue

            # DOC: Get the alert method within the water depth area
            wd_alert_method = gdf_alert_method_3857.cx[wd_bounds[0]:wd_bounds[2], wd_bounds[1]:wd_bounds[3]]    # !!!: There is no an explicit relation between alert area and method, so we look at which method is contained in the alert area
            if wd_alert_method.empty:
                continue    # DOC: this should only when area is related to a not-flooding alert method or if method is not acoustic

            wd_alert_method = wd_alert_method.drop_duplicates(subset='id')
            
            # DOC: Collect the alert method GeoDataFrame
            wd_alert_method = wd_alert_method.to_crs(kwargs['t_srs'])
            alert_method_gdfs.append(wd_alert_method)
            


            # DOC: Mark (flooded) buildings in the alert area with relative alert method
            alert_buildings = gdf_buildings_3857[gdf_buildings_3857.geometry.intersects(wd_alert_area.geometry.union_all())]    
            # ???: Filter only alert_buildings.is_flooded → (logically yes but maybe it is important to have the alert method for all buildings in the alert area)
            if alert_buildings.empty:
                continue
            wd_alert_method['location'] = wd_alert_method['geometry'].apply(lambda geom: [geom.x, geom.y])
            wd_alert_method = pd.DataFrame(wd_alert_method.drop(columns='geometry'))
            wd_alert_method = wd_alert_method[self._alert_method_basic_attributes + ['location']]
            wd_alert_method = _utils.df_dt_col_to_isoformat(wd_alert_method)
            wd_alert_method = json.loads(wd_alert_method.to_json(orient='records'))
            for i_building, _ in alert_buildings.iterrows():
                gdf_buildings.at[i_building, 'alert_method'] = wd_alert_method


        if alert_method_gdfs:
            alert_method_gdf = pd.concat(alert_method_gdfs, ignore_index=True).drop_duplicates(subset='id').reset_index(drop=True)
            alert_method_gdf = _utils.df_dt_col_to_isoformat(alert_method_gdf)
            alert_method_collection = alert_method_gdf.to_geo_dict()
        else:
            alert_method_collection = { 'type': 'FeatureCollection', 'features': [] }

        # DOC: Return the updated building GeoDataFrame with the alert method and the alert method collection    
        return gdf_buildings, alert_method_collection





_ADD_OPS = {

    # DOC: Additional operations for VENEZIA_WFS_PROVIDER
    _consts._VENEZIA_WFS_PROVIDER: {
        NearbyPumps.name: NearbyPumps,
        AlertMethod.name: AlertMethod
    }

}


def get_ops_by_provider(provider: str) -> dict[str | AdditionalOperation]:
    """
    Retrieve the additional operation by provider and operation name.
    """
    
    if provider in _ADD_OPS:
        return _ADD_OPS[provider]

    if provider.startswith(_consts._RER_REST_PROVIDER):
        provider = _consts._RER_REST_PROVIDER
    elif provider.startswith(_consts._VENEZIA_WFS_PROVIDER):
        provider = _consts._VENEZIA_WFS_PROVIDER
        
    return _ADD_OPS.get(provider, [])


def get_op_by_name(provider: str, op_name: str) -> AdditionalOperation | None:
    """
    Retrieve the additional operation by provider and operation name.
    """
    
    ops = get_ops_by_provider(provider)
    if op_name in ops:
        return ops[op_name]
    
    return None