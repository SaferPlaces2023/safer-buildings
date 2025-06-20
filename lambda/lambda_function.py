from safer_buildings import parse_event
from safer_buildings import compute_flood as main_function


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
