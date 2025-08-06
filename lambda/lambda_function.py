import pprint
from safer_buildings import _consts
from safer_buildings import parse_event
from safer_buildings import compute_flood as main_function


# DOC: An example test input for directed project area:
# {
#   "water": "s3://saferplaces.co/Directed/data-fabric-rwl2/WD_radar_saferplaces_edekmzo.tif",
#   "building": "s3://saferplaces.co/Directed/process_out/SaferBuildingsService/Data/buildings-default-area__rer-rest_overture.geojson",
#   "out": "s3://saferplaces.co/Directed/process_out/SaferBuildingsService/test-rimini-wd-00-output-00.geojson",
#   "provider": "RER-REST/1",
#   "summary": true
# }


def lambda_handler(event, context):
    """
    lambda_handler - lambda function
    """
    kwargs = parse_event(event, main_function)

    result = main_function(**kwargs)
    status_code = result.pop(_consts._ERROR_CODE_KEY, _consts._OUTPUT_SUCCESS_CODE)

    return {
        "statusCode": status_code, 
        "body": {
            "result": result   
        }
    }


if __name__ == "__main__":
    # event = {
    #     "water": "s3://saferplaces.co/packages/safer_rain/CLSA_LiDAR/CLSA_LiDAR.tif",
    #     "buildings": "",
    #     "debug": "true"
    # }

    # kwargs = parse_event(event, main_function)
    # main_function(**kwargs)

    event = {
        "water": "INVALIDO!", #"s3://saferplaces.co/Venezia/WaterDepthsv2/ICON_2I_SURFACE_PRESSURE_LEVELS_tp/2025-07-28/00-00/water_depth_bacino2_forecast_acc_12h_2025-07-28_00-00_01h-12h.tif", #"s3://saferplaces.co/Safer-Buildings/test/venezia-wd-400mm-1h.tif"
        "building": "s3://saferplaces.co/Venezia/shapes/buildings/building_2.shp", #,
        "wd_thresh": None,
        "bbox": None,
        "out": "s3://saferplaces.co/Venezia/SaferBuildings/water_depth_bacino2_forecast_acc_12h_2025-07-28_00-00_01h-12h__building_2_nearby-pumps.geojson",  #"s3://saferplaces.co/Safer-Buildings/test/venezia-wd-400mm-1h-flood-buildings-add-ops.geojson",
        "t_srs": "EPSG:4326",
        "provider": f'{_consts._VENEZIA_WFS_PROVIDER}', #/v_pc_p0106011_scuole',
        "filters": None,
        "only_flood": False,
        "stats": False,
        "summary": False,
        "summary_on": "subtype",    # None,
        "add_ops": None, #{
        #     module_add_ops.NearbyPumps.name: {
        #         "wd_buffer": 0.0,
        #     },
        #     module_add_ops.AlertMethod.name: {
        #         "wd_buffer": 0.0,
        #     }
        # },
        "out_geojson": False,

        "version": False,
        "debug": True,
        "verbose": False,
    }

    output_lambda = lambda_handler(event, None)

    print('\n', '=' * 50, '\n')
    print("Output from lambda_handler", "\n")
    pprint.pprint(output_lambda)
    print('\n', '=' * 50, '\n')
