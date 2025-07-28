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

    name = 'nearby_pumps'
    args = [
        # DOC: Maximum distance to consider a pump as "nearby" (in meters). Default is 1000.0 meters.
        'max_distance',
        # DOC: If True, only consider pumps that intersect with the water depth area of the building. Default is False.
        'wd_intersect'              
    ]

    _pumps_layer_id = 'mv_risorse_p0109103_pompe'
    _nearby_pumps_basic_attributes = ['id', 'modello', 'distance', 'location']

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


        nearby_pumps_gdfs = []
        # DOC: Foreach building, retrieve (if any) the pump within the water depth area
        gdf_buildings['nearby_pumps'] = [dict() for _ in range(len(gdf_buildings))]
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
            
            # DOC: Add pump information to the building
            nearby_pumps = nearby_pumps.to_crs(kwargs['t_srs'])
            nearby_pumps['location'] = nearby_pumps['geometry'].apply(lambda geom: [geom.x, geom.y])
            
            nearby_pumps_gdfs.append(nearby_pumps)
            
            nearby_pumps = pd.DataFrame(nearby_pumps.drop(columns='geometry'))
            nearby_pumps = nearby_pumps[self._nearby_pumps_basic_attributes]
            nearby_pumps = _utils.df_dt_col_to_isoformat(nearby_pumps)
            nearby_pumps = json.loads(nearby_pumps.to_json(orient='records'))
            gdf_buildings.at[i_building, 'nearby_pumps'] = nearby_pumps

        # DOC: Return the updated GeoDataFrame with the nearest pumps

        if nearby_pumps_gdfs:
            nearby_pumps_gdf = pd.concat(nearby_pumps_gdfs, ignore_index=True).drop_duplicates(subset='id').reset_index(drop=True)
            nearby_pumps_gdf = _utils.df_dt_col_to_isoformat(nearby_pumps_gdf)
            nearby_pumps_collection = nearby_pumps_gdf.to_geo_dict()
        else:
            nearby_pumps_collection = { 'type': 'FeatureCollection', 'features': [] }

        return gdf_buildings, nearby_pumps_collection


        
_ADDITIONAL_OPERATIONS = {

    # DOC: Additional operations for VENEZIA_WFS_PROVIDER
    _consts._VENEZIA_WFS_PROVIDER: {

        NearbyPumps.name: NearbyPumps

    }

}


def get_ops_by_provider(provider: str) -> dict[str | AdditionalOperation]:
    """
    Retrieve the additional operation by provider and operation name.
    """
    
    if provider in _ADDITIONAL_OPERATIONS:
        return _ADDITIONAL_OPERATIONS[provider]

    if provider.startswith(_consts._RER_REST_PROVIDER):
        provider = _consts._RER_REST_PROVIDER
    elif provider.startswith(_consts._VENEZIA_WFS_PROVIDER):
        provider = _consts._VENEZIA_WFS_PROVIDER
        
    return _ADDITIONAL_OPERATIONS.get(provider, [])


def get_op_by_name(provider: str, op_name: str) -> AdditionalOperation | None:
    """
    Retrieve the additional operation by provider and operation name.
    """
    
    ops = get_ops_by_provider(provider)
    if op_name in ops:
        return ops[op_name]
    
    return None