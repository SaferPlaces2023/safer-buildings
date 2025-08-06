import os
import unittest

from safer_buildings import compute_flood as main_python
from safer_buildings import module_s3, filesystem, _consts


class Test(unittest.TestCase):
    """
    Test class for the greenamptr module.
    """

    def test_rimini_001(self):
        f"""
        CLI Command:
        safer-buildings --water s3://saferplaces.co/Safer-Buildings/test/rimini-wd-100mm-1h.tif --out s3://saferplaces.co/Safer-Buildings/test/rimini-wd-100mm-1h-flood-buildings.geojson --provider {_consts._RER_REST_PROVIDER} --summary --stats --debug
        """
        
        args = {
            "water": "s3://saferplaces.co/Safer-Buildings/test/rimini-wd-100mm-1h.tif",
            "building": None,
            "wd_thresh": None,
            "bbox": None,
            "out": "s3://saferplaces.co/Safer-Buildings/test/rimini-wd-100mm-1h-flood-buildings.geojson",
            "t_srs": None,
            "provider": _consts._RER_REST_PROVIDER,
            "filters": None,
            "only_flood": False,
            "stats": True,
            "summary": True,
            "summary_on": None,
            "out_geojson": False,

            "version": False,
            "debug": True,
            "verbose": False,
        }

    
        result = main_python( ** args )

        try:
            self.assertEqual(result.pop(_consts._ERROR_CODE_KEY, _consts._OUTPUT_SUCCESS_CODE), _consts._OUTPUT_SUCCESS_CODE)
        except:
            print(f"Test 'test_vimercate_001' failed. Error: {result}")

        fileout = module_s3.s3_download(
            uri = args['out'],
            fileout = filesystem.tempfilename(suffix='.geojson', prefix='test-rimini-001-')
        )

        try:
            self.assertTrue(os.path.exists(fileout))
            print(f"Test 'test_rimini_001' successfully completed.")
        except:
            print(f"Test 'test_rimini_001' failed. Output file not found: {fileout}")
        

if __name__ == '__main__':
    unittest.main()
