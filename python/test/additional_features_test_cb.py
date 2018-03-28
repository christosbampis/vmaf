__copyright__ = "Copyright 2016-2018, Netflix, Inc."
__license__ = "Apache, Version 2.0"

import os
import unittest
import re

from vmaf.config import VmafConfig
from vmaf.core.feature_extractor import VmafFeatureExtractor, SpEEDMatlabFeatureExtractor, VmafDiffFeatureExtractor
from vmaf.core.asset import Asset


class FeatureExtractorTest(unittest.TestCase):

    def tearDown(self):
        if hasattr(self, 'fextractor'):
            self.fextractor.remove_results()
        pass

    def test_run_vmaf_fextractor(self):
        print 'test on running VMAF feature extractor...'
        ref_path = VmafConfig.test_resource_path("yuv", "src01_hrc00_576x324.yuv")
        dis_path = VmafConfig.test_resource_path("yuv", "src01_hrc01_576x324.yuv")
        asset = Asset(dataset="test", content_id=0, asset_id=0,
                      workdir_root=VmafConfig.workdir_path(),
                      ref_path=ref_path,
                      dis_path=dis_path,
                      asset_dict={'width':576, 'height':324})

        asset_original = Asset(dataset="test", content_id=0, asset_id=1,
                      workdir_root=VmafConfig.workdir_path(),
                      ref_path=ref_path,
                      dis_path=ref_path,
                      asset_dict={'width':576, 'height':324})

        self.fextractor = VmafFeatureExtractor(
            [asset, asset_original],
            None, fifo_mode=True,
            result_store=None
        )
        self.fextractor.run()

        results = self.fextractor.results

        self.assertAlmostEqual(results[0]['VMAF_feature_vif_score'], 0.4460930625, places=4)
        self.assertAlmostEqual(results[0]['VMAF_feature_motion_score'], 4.04982535417, places=4)
        self.assertAlmostEqual(results[0]['VMAF_feature_motion2_score'], 3.8953518541666665, places=4)
        self.assertAlmostEqual(results[0]['VMAF_feature_adm_score'],0.93458780728708746, places=4)
        self.assertAlmostEqual(results[0]['VMAF_feature_adm2_score'],0.93458780728708746, places=4) # at version 0.2.4b (ioannis adm fix), adm and adm2 are now identical
        self.assertAlmostEqual(results[0]['VMAF_feature_ansnr_score'], 23.5095715208, places=4)

        self.assertAlmostEqual(results[0]['VMAF_feature_vif_num_score'], 712650.023478, places=0)
        self.assertAlmostEqual(results[0]['VMAF_feature_vif_den_score'], 1597314.95249, places=0)
        self.assertAlmostEqual(results[0]['VMAF_feature_adm_num_score'], 371.83541406249998, places=4)
        self.assertAlmostEqual(results[0]['VMAF_feature_adm_den_score'], 397.83378972916671, places=4)
        self.assertAlmostEqual(results[0]['VMAF_feature_anpsnr_score'], 34.164776875, places=4)

        self.assertAlmostEqual(results[0]['VMAF_feature_vif_scale0_score'], 0.363420489439, places=4)
        self.assertAlmostEqual(results[0]['VMAF_feature_vif_scale1_score'], 0.766647542135, places=4)
        self.assertAlmostEqual(results[0]['VMAF_feature_vif_scale2_score'], 0.862854666902, places=4)
        self.assertAlmostEqual(results[0]['VMAF_feature_vif_scale3_score'], 0.915971778036, places=4)

        self.assertAlmostEqual(results[0]['VMAF_feature_adm_scale0_score'], 0.90791933424090698, places=4)
        self.assertAlmostEqual(results[0]['VMAF_feature_adm_scale1_score'], 0.89395660507453423, places=4)
        self.assertAlmostEqual(results[0]['VMAF_feature_adm_scale2_score'], 0.93010045185439161, places=4)
        self.assertAlmostEqual(results[0]['VMAF_feature_adm_scale3_score'], 0.96503534602850938, places=4)

        self.assertAlmostEqual(results[0]['VMAF_feature_vif2_score'], 0.72722361912801026, places=4)
        self.assertAlmostEqual(results[0]['VMAF_feature_adm3_score'], 0.92425293429958566, places=4)

        self.assertAlmostEqual(results[1]['VMAF_feature_vif_score'], 1.0, places=4)
        self.assertAlmostEqual(results[1]['VMAF_feature_motion_score'], 4.04982535417, places=4)
        self.assertAlmostEqual(results[1]['VMAF_feature_motion2_score'], 3.8953518541666665, places=4)
        self.assertAlmostEqual(results[1]['VMAF_feature_adm_score'], 1.0, places=4)
        self.assertAlmostEqual(results[1]['VMAF_feature_adm2_score'], 1.0, places=4)
        self.assertAlmostEqual(results[1]['VMAF_feature_ansnr_score'], 31.2714392708, places=4)

        self.assertAlmostEqual(results[1]['VMAF_feature_vif_num_score'], 1597314.86733, places=0)
        self.assertAlmostEqual(results[1]['VMAF_feature_vif_den_score'], 1597314.95249, places=0)
        self.assertAlmostEqual(results[1]['VMAF_feature_adm_num_score'], 397.83378972916671, places=4)
        self.assertAlmostEqual(results[1]['VMAF_feature_adm_den_score'], 397.83378972916671, places=4)
        self.assertAlmostEqual(results[1]['VMAF_feature_anpsnr_score'], 41.9266444375, places=4)

        self.assertAlmostEqual(results[1]['VMAF_feature_vif_scale0_score'], 1.0, places=4)
        self.assertAlmostEqual(results[1]['VMAF_feature_vif_scale1_score'], 1.0, places=4)
        self.assertAlmostEqual(results[1]['VMAF_feature_vif_scale2_score'], 1.0, places=4)
        self.assertAlmostEqual(results[1]['VMAF_feature_vif_scale3_score'], 1.0, places=4)

        self.assertAlmostEqual(results[1]['VMAF_feature_adm_scale0_score'], 1.0, places=4)
        self.assertAlmostEqual(results[1]['VMAF_feature_adm_scale1_score'], 1.0, places=4)
        self.assertAlmostEqual(results[1]['VMAF_feature_adm_scale2_score'], 1.0, places=4)
        self.assertAlmostEqual(results[1]['VMAF_feature_adm_scale3_score'], 1.0, places=4)

        self.assertAlmostEqual(results[1]['VMAF_feature_vif2_score'], 1.0, places=4)
        self.assertAlmostEqual(results[1]['VMAF_feature_adm3_score'], 1.0, places=4)

    def test_run_speed_matlab_fextractor(self):

        print 'test on running SpEED (Matlab) feature extractor...'

        ref_path = VmafConfig.test_resource_path("yuv", "src01_hrc00_576x324.yuv")
        dis_path = VmafConfig.test_resource_path("yuv", "src01_hrc01_576x324.yuv")

        asset = Asset(dataset="test", content_id=0, asset_id=0,
                      workdir_root=VmafConfig.workdir_path(),
                      ref_path=ref_path,
                      dis_path=dis_path,
                      asset_dict={'width': 576, 'height': 324})

        asset_original = Asset(dataset="test", content_id=0, asset_id=1,
                      workdir_root=VmafConfig.workdir_path(),
                      ref_path=ref_path,
                      dis_path=ref_path,
                      asset_dict={'width': 576, 'height': 324})

        self.fextractor = SpEEDMatlabFeatureExtractor(
            [asset, asset_original],
            None, fifo_mode=True,
            result_store=None
        )
        self.fextractor.run(parallelize=False)

        results = self.fextractor.results

        # S-SpEED assertions on first frame

        self.assertAlmostEqual(results[0].result_dict[self.fextractor.TYPE + '_sspeed_2_scores'][0], 13.510418, places=4)
        self.assertAlmostEqual(results[0].result_dict[self.fextractor.TYPE + '_sspeed_3_scores'][0], 7.211881, places=4)
        self.assertAlmostEqual(results[0].result_dict[self.fextractor.TYPE + '_sspeed_4_scores'][0], 4.921501, places=4)

        # T-SpEED assertions on third frame

        self.assertAlmostEqual(results[0].result_dict[self.fextractor.TYPE + '_tspeed_2_scores'][2], 32.994605, places=4)
        self.assertAlmostEqual(results[0].result_dict[self.fextractor.TYPE + '_tspeed_3_scores'][2], 22.404285, places=4)
        self.assertAlmostEqual(results[0].result_dict[self.fextractor.TYPE + '_tspeed_4_scores'][2], 15.233468, places=4)

    def test_run_vmaf_diff_fextractor(self):
        print 'test on running VMAF diff feature extractor...'
        ref_path = VmafConfig.test_resource_path("yuv", "src01_hrc00_576x324.yuv")
        dis_path = VmafConfig.test_resource_path("yuv", "src01_hrc01_576x324.yuv")
        asset = Asset(dataset="test", content_id=0, asset_id=0,
                      workdir_root=VmafConfig.workdir_path(),
                      ref_path=ref_path,
                      dis_path=dis_path,
                      asset_dict={'width':576, 'height':324})

        asset_original = Asset(dataset="test", content_id=0, asset_id=1,
                      workdir_root=VmafConfig.workdir_path(),
                      ref_path=ref_path,
                      dis_path=ref_path,
                      asset_dict={'width':576, 'height':324})

        self.fextractor = VmafDiffFeatureExtractor(
            [asset, asset_original],
            None, fifo_mode=True,
            result_store=None
        )
        self.fextractor.run()

        results = self.fextractor.results

        # T-VIF assertions on THIRD frame

        self.assertAlmostEqual(results[0]['VMAFDiff_feature_vif_diff_scores'][2], 0.287009, places=4)
        self.assertAlmostEqual(results[0]['VMAFDiff_feature_vif_scale0_diff_scores'][2], 0.233937484911, places=4)


if __name__ == '__main__':
    unittest.main()
