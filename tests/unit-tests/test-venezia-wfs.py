import os
import unittest

from safer_buildings import compute_flood as main_python, module_add_ops
from safer_buildings import module_s3, filesystem, _consts


class Test(unittest.TestCase):
    """
    Test class for the greenamptr module.
    """

    # TODO: Static files to be added in handling of Venezia WFS provider
    # - DataFrame of WFS Services
        # - Full builings: Overtoure OPPURE
        # - Full buildings Venezia WFS provider (subset of the layers , ask fischietti which are the ones of interest) (in check-args crop to bbox)
    # - Full Pumps gdf (in add ops-crop to bbox)
    # - Full Alert area gdf (in add ops-crop to bbox)
    # - Full Alert method gdf (in add ops-crop to bbox)

    def test_venezia_001(self):
        f"""
        CLI Command:
        safer-buildings --water s3://saferplaces.co/Safer-Buildings/test/venezia-wd-400mm-1h.tif --out s3://saferplaces.co/Safer-Buildings/test/venezia-wd-400mm-1h-flood-buildings-nearest_pump.geojson --provider VENEZIA-WFS/v_pc_p0106011_scuole --summary --stats --add_ops nearby_pumps --debug
        """
        
        bbox_lotto_1 = [11.9360964, 44.9937129, 12.3501891, 45.3141472]

        args_big = {
            "water": "s3://saferplaces.co/Venezia/WaterDepthsv2/ICON_2I_SURFACE_PRESSURE_LEVELS_tp/2025-07-28/00-00/water_depth_bacino2_forecast_acc_12h_2025-07-28_00-00_01h-12h.tif", #"s3://saferplaces.co/Safer-Buildings/test/venezia-wd-400mm-1h.tif"
            # "building": "s3://saferplaces.co/Venezia/shapes/buildings/building_2.shp", #,
            "wd_thresh": None,
            "bbox": bbox_lotto_1,
            "out": "s3://saferplaces.co/Venezia/SaferBuildings/geojson/water_depth_bacino2_forecast_acc_12h_2025-07-28_00-00_01h-12h__building_2_outformat.geojson",  #"s3://saferplaces.co/Safer-Buildings/test/venezia-wd-400mm-1h-flood-buildings-add-ops.geojson",
            "t_srs": "EPSG:4326",
            # "provider": f'{_consts._VENEZIA_WFS_PROVIDER}', #/v_pc_p0106011_scuole',
            "provider": f'{_consts._VENEZIA_WFS_CRITICAL_SITES_PROVIDER}',
            "filters": None,
            "only_flood": False,
            "stats": False,
            "summary": False,
            "summary_on": "subtype",    # None,
            "add_ops": {
                module_add_ops.NearbyPumps.name: {
                    # "wd_buffer": 2000.0,
                    "max_distance": 2000.0,  # DOC: Set the maximum distance for nearby pumps
                },
                module_add_ops.AlertMethod.name: {
                    "wd_buffer": 200.0,
                }
            },
            "out_geojson": False,

            "version": False,
            "debug": True,
            "verbose": False,
        }

        args_small = {
            "water":"s3://saferplaces.co/Safer-Buildings/test/venezia-wd-400mm-1h.tif",
            "building": None, #"s3://saferplaces.co/Venezia/shapes/buildings/building_2.shp", #,
            "wd_thresh": None,
            "bbox": None,
            "out": "s3://saferplaces.co/Safer-Buildings/test/venezia-wd-400mm-1h-flood-buildings-flooded-prj.geojson",
            "t_srs": "EPSG:4326",
            "provider": f'{_consts._VENEZIA_WFS_CRITICAL_SITES_PROVIDER}',
            "filters": None,
            "only_flood": False,
            "stats": False,
            "summary": False,
            "summary_on": "subtype",    # None,
            "add_ops": {
                module_add_ops.NearbyPumps.name: {
                    "max_distance": 4000.0,
                },
                module_add_ops.AlertMethod.name: {
                    "wd_buffer": 2000.0,
                }
            },
            "out_geojson": False,

            "version": False,
            "debug": True,
            "verbose": False,
        }

    
        args = args_small
        
        result = main_python( ** args )

        try:
            self.assertEqual(result.pop(_consts._ERROR_CODE_KEY, _consts._OUTPUT_SUCCESS_CODE), _consts._OUTPUT_SUCCESS_CODE)
        except:
            print(f"Test 'test_vimercate_001' failed. Error: {result}")
        
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
