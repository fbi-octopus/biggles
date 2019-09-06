import sys
import unittest

import biggles.simulate as sim
import numpy as np

class TSTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def test_control_parameters(self):
        self.assertTrue("parameters" in sim.SIMULATION)
        self.assertTrue("original" in sim.SIMULATION)
        self.assertTrue("options" in sim.SIMULATION)
        self.assertTrue("two_mode_tracks" in sim.SIMULATION)
        self.assertTrue("constraining_regions" in sim.SIMULATION)
        
    def test_obs_cov(self):
        sim.SIMULATION['parameters']['observation_probability'] = 1.0
        sim.SIMULATION['parameters']['observation_covariance'] = [
            [0.001, 0.0],[0.0, 0.001]
        ]
        sim.SIMULATION['parameters']['survival_probability'] = .99
        sim.SIMULATION['options']['creation_events'] = 'uniform'
        sim.SIMULATION['options']['split_probability'] = 0.0
        sim.SIMULATION['options']['zig_zag'] = 0.0
        sim.SIMULATION['options']['directionality_proportion'] = 0.0
        sim.SIMULATION['options']['position_sigma'] = 0.0
        gt = sim.generate_tracks()
        tracks, _, _ = sim.observation_model(gt)
        cms = []
        for track in tracks:
            if len(track['observations'])>1:
                obs = np.array(track['observations'])[:,0:2]
                cms.append(np.cov(obs.T))
        cm = np.mean(cms, axis = 0)
        self.assertAlmostEqual(cm[0,0], 0.001, delta = 0.0001)
        self.assertAlmostEqual(cm[1,1], 0.001, delta = 0.0001)
        self.assertEqual(cm[0,1], cm[1,0])
        self.assertLess(cm[0,1], 0.0001)
        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main(argv=sys.argv[:1])
