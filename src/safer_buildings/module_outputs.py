import os
import json

import geopandas as gpd

from . import _consts, _utils, filesystem, module_s3
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
        'buildings_count': len(flooded_buildings),
        'flooded_buildings_count': int(flooded_buildings[_consts._COL_IS_FLOODED].sum()),
        **( {'provider': provider} if provider else dict() ),
        ** summary_stats,
        ** add_ops_output_data    # !!!: Fischietti ci dirà se tenere un'intera feature collection nei metadati, in mia opinione eccessivo, meglio avere più file separati <out-name>.<add-op-name>.<out-format> (comunque già implementato)
    }
    feature_collection = _utils.set_crs_feature_collection(feature_collection, t_srs)
    
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
    
    # DOC: support variable used to handle both local and S3 outputs
    local_output = out
    # DOC: Get ext to determine the output file format, its drivers and auxiliary files (if any)
    out_ext = filesystem.justext(out)

    # DOC: Postifx is used for addtitional files, e.g. additional operations output
    if out_postfix:
        out = filesystem.forceext(out, f"{out_postfix}.{out_ext}")

    # DOC: Read the feature collection and convert to GeoDataFrame (useful for saving file in different formats)
    gdf_fc = gpd.GeoDataFrame.from_features(feature_collection['features'], crs=feature_collection['crs']['properties']['name'])

    # DOC: If the output is S3, local output is a temporary file
    if out.startswith('s3://'):
        local_output = _utils.temp_filename(ext=out_ext, prefix='safer-buildings_out')
    
    # DOC: Save the GeoDataFrame to the local output file
    gdf_fc.to_file(filename=local_output, driver=filesystem._GPD_DRIVERS(out_ext))
    
    # DOC: If the output is GeoJSON, restore its metadata
    if out_ext == 'geojson' and 'metadata' in feature_collection:
        with open(local_output, 'r', encoding='utf8') as f:
            data = json.load(f)
            data['metadata'] = feature_collection['metadata']
        with open(local_output, 'w', encoding='utf8') as f:
            json.dump(data, f, indent=2)
    
    # DOC: Keep track of auxiliary files based on the output format
    local_output = [local_output] + filesystem.get_aux_files(local_output)

    # DOC: If the output is S3, upload the local output file(s) to S3
    if out.startswith('s3://'):
        for lout in local_output:
            module_s3.s3_upload(filename=lout, uri=filesystem.forceext(out, filesystem.justext(lout)))

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