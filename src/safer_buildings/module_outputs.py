import os
import json

import geopandas as gpd

from . import _utils, filesystem, module_s3
from .module_log import Logger



def prepare_feture_collection(
    flooded_buildings: gpd.GeoDataFrame,
    t_srs: str,
    provider: str,
    summary_stats: dict = dict(),
    add_ops_output_data: dict = dict()
):
    """
    Prepare a GeoJSON FeatureCollection from the flooded buildings data.
    """
    
    flooded_buildings = flooded_buildings.to_crs(t_srs)
    flooded_buildings = _utils.safe_json_df(flooded_buildings)
    feature_collection = flooded_buildings.to_geo_dict()
    
    feature_collection['metadata'] = {
        'provider': provider,
        'buildings_count': len(flooded_buildings),
        'flooded_buildings_count': int(flooded_buildings['is_flooded'].sum()),
        ** summary_stats,
        ** add_ops_output_data
    }
    feature_collection['crs'] = {
        "type": "name",
        "properties": {
            "name": f"urn:ogc:def:crs:{t_srs.replace(':', '::')}"  # REF: https://gist.github.com/sgillies/1233327 lines 256:271
        }
    }
    
    Logger.debug(f"## Buildings feature collection prepared with {len(feature_collection['features']) if 'features' in feature_collection else 0} features.")
    
    return feature_collection


def save_results(
    feature_collection: dict,
    out: str,
    out_postfix: str = None
):
    """
    Save the results to a file or upload to S3.
    """
    
    if out_postfix:
        out_ext = filesystem.justext(out)
        out = out.replace(out_ext, f"{out_postfix}.{out_ext}")

    if out.startswith('s3://'):
        out_tmp = _utils.temp_filename(ext='geojson', prefix='safer-buildings_out')
        with open(out_tmp, 'w') as f:
            json.dump(feature_collection, f, indent=2)
        module_s3.s3_upload(filename = out_tmp, uri = out)    
    else:
        with open(out, 'w') as f:
            json.dump(feature_collection, f, indent=2)

    Logger.debug(f"## Results saved to {out}")


def prepare_output(
    feature_collection: dict,
    out: str,
    out_geojson: bool,
    compute_summary: bool
):
    output = feature_collection
    if not out_geojson:
        geojson_ref_key = 's3_uri' if out.startswith('s3://') else 'geojson_file'
        summary_info = {'summary': output['metadata']['summary']} if compute_summary else dict()
        output = {
            geojson_ref_key: out,
            ** summary_info
        }

    Logger.debug(f"## Output prepared. Main fields are {list(output.keys())}.")
    
    return output