import os
import unittest

from safer_buildings import compute_flood as main_python
from safer_buildings import module_s3


class Test(unittest.TestCase):
    """
    Test class for the greenamptr module.
    """

    def test_rimini_00(self):
        """
        CLI Command:
        safer-buildings --water s3://saferplaces.co/Directed/process_out/SaferBuildingsService/rimini-wd.tif --building s3://saferplaces.co/Directed/process_out/SaferBuildingsService/Data/buildings-default-area__rer-rest_overture.geojson --out s3://saferplaces.co/Directed/process_out/SaferBuildingsService/rimini-wd-buildings.geojson --provider RER-REST/1 --summary
        """
        
        args = {
            'water':        "s3://saferplaces.co/Directed/process_out/SaferBuildingsService/rimini-wd.tif",
            'building':     "s3://saferplaces.co/Directed/process_out/SaferBuildingsService/Data/buildings-default-area__rer-rest_overture.geojson",
            'out':          "s3://saferplaces.co/Directed/process_out/SaferBuildingsService/test-rimini-wd-00-output-00.geojson",
            'provider':     "RER-REST/1",
            'summary':      True,
            
            # 'debug':        False,
            # 'verbose':      True,
        }

    
        main_python( ** args )
        
        
        fileout = module_s3.s3_download(
            uri = args['out'],
            fileout = 'test-rimini-01-output.geojson',
        )

        self.assertTrue(os.path.exists(fileout))    # TODO: We should check its content.
        

if __name__ == '__main__':
    unittest.main()
