import os
import unittest

from safer_buildings import compute_flood as main_python
from safer_buildings import module_s3, filesystem, _consts


class Test(unittest.TestCase):
    """
    Test class for the greenamptr module.
    """

    def test_venezia_001(self):
        f"""
        CLI Command:
        safer-buildings --water s3://saferplaces.co/Safer-Buildings/test/venezia-wd-400mm-1h.tif --out s3://saferplaces.co/Safer-Buildings/test/venezia-wd-400mm-1h-flood-buildings.geojson --provider {_consts._VENEZIA_WFS_PROVIDER} --summary --stats --debug
        """
        
        args = {
            "water": "s3://saferplaces.co/Safer-Buildings/test/venezia-wd-400mm-1h.tif",
            "building": None,
            "wd_thresh": None,
            "bbox": None,
            "out": "s3://saferplaces.co/Safer-Buildings/test/venezia-wd-400mm-1h-flood-buildings.geojson",
            "t_srs": None,
            "provider": _consts._VENEZIA_WFS_PROVIDER,
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

    
        main_python( ** args )
        
        fileout = module_s3.s3_download(
            uri = args['out'],
            fileout = filesystem.tempfilename(suffix='.geojson', prefix='test-venezia-001-')
        )

        try:
            self.assertTrue(os.path.exists(fileout))
            print(f"Test 'test_venezia_001' successfully completed.")
        except:
            print(f"Test 'test_venezia_001' failed. Output file not found: {fileout}")
        

if __name__ == '__main__':
    unittest.main()
