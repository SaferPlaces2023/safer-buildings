from safer_buildings import parse_event
from safer_buildings import compute_flood, main as main_function


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

    res = main_function(**kwargs)

    return {
        "statusCode": 200, 
        "body": {
            "result": res   
        }
    }


# if __name__ == "__main__":
#     event = {
#         "water": "s3://saferplaces.co/packages/safer_rain/CLSA_LiDAR/CLSA_LiDAR.tif",
#         "buildings": "",
#         "debug": "true"
#     }

#     kwargs = parse_event(event, main_function)
#     main_function(**kwargs)
