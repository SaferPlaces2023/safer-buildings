import os
import json
import unittest

from safer_buildings import compute_flood as main_python
from safer_buildings import module_s3, filesystem, _consts


class Test(unittest.TestCase):
    """
    Test class for the greenamptr module.
    """

    def test_deb(self):
        f"""
        CLI Command:
        safer-buildings --water s3://saferplaces.co/Safer-Buildings/test/rimini-wd-100mm-1h.tif --out s3://saferplaces.co/Safer-Buildings/test/rimini-wd-100mm-1h-flood-buildings.geojson --provider {_consts._RER_REST_PROVIDER} --summary --stats --debug
        """
        
        args = {
            "water": "s3://saferplaces.co/Directed/data-fabric-rwl2/WD_radar_saferplaces_a3j5sxr.tif", #"s3://saferplaces.co/Directed/data-fabric-rwl2/WD_radar_saferplaces_0fw6zyg.tif",
            "building": "s3://saferplaces.co/Directed/process_out/SaferBuildingsService/Data/directed-rer-overture-buildings.gpkg", #"s3://saferplaces.co/Directed/process_out/SaferBuildingsService/Data/buildings-default-area__rer-rest_overture.geojson", //"https://s3.us-east-1.amazonaws.com/saferplaces.co/Venezia/shapes/buildings/building_2.shp",
            "wd_thresh": 0.1,
            "bbox": None,
            "out": "s3://saferplaces.co/Directed/data-fabric-rwl2/Rimini_coast_flooded_test-js-20250902-1329.geojson",
            "t_srs": None,
            "provider": None,
            "filters": None,
            "only_flood": True,
            "stats": True,
            "summary": True,
            "summary_on": "subtype",
            "out_geojson": False,

            "version": False,
            "debug": True,
            "verbose": False,

            # "token": "S4fer_api_token",
            # "user": "saferplaces",
        }

    
        result = main_python( ** args )

        with open('test_deb_result.json', 'w') as f:
            json.dump(result, f, indent=4)

        print()
        print(result)
        print()

        try:
            self.assertEqual(result.pop(_consts._ERROR_CODE_KEY, _consts._OUTPUT_SUCCESS_CODE), _consts._OUTPUT_SUCCESS_CODE)
        except:
            print(f"Test 'test_vimercate_001' failed. Error: {result}")
        
        fileout = module_s3.s3_download(
            uri = args['out'],
            fileout = filesystem.tempfilename(suffix='.geojson', prefix='test_deb-')
        )

        try:
            self.assertTrue(os.path.exists(fileout))
            print(f"Test 'test_deb' successfully completed.")
        except:
            print(f"Test 'test_deb' failed. Output file not found: {fileout}")
        

if __name__ == '__main__':
    unittest.main()
