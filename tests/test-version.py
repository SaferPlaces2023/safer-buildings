import os
import unittest

from safer_buildings import compute_flood as main_python


class Test(unittest.TestCase):
    """
    Test class for the greenamptr module.
    """

    def test_version(self):
        """
        CLI Command:
        safer-buildings --version
        """
        
        args = {
            'version': True
        }

    
        main_python( ** args )

        self.assertTrue(True)  # Just checking if the version command runs without errors.
        

if __name__ == '__main__':
    unittest.main()
