import unittest
import os
from safer_buildings import compute_flood as main_python

class Test(unittest.TestCase):
    """
    Test class for the greenamptr module.
    """

    def test_CLSA_LiDAR(self):
        """
        test
        """
        water = "s3://saferplaces.co/packages/safer_buildings/CLSA_LiDAR/CLSA_LiDAR.wd.tif"

    
        main_python(...
                        debug=False, verbose=True)

        self.assertTrue(os.path.exists(fileout))
        
      


if __name__ == '__main__':
    unittest.main()
