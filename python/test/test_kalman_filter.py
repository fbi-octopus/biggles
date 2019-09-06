import json
import StringIO
import os
import unittest
import sys

from biggles import ModelParameters, Partition, ObservationCollection, TrackCollection

test_data_path = sys.argv[1]

class PartitionTestCase(unittest.TestCase):
    def setUp(self):
        self.p = Partition()

        p_json = open(os.path.join(test_data_path, 'simple-ground-truth.json')).read()
        self.simple = Partition.from_json(p_json)

        self.mp = ModelParameters.from_json(json.dumps({
            'births_per_frame': 0.0,
            'clutter_per_frame': 0.0,
            'survival_probability': 0.0,
            'observation_probability': 0.0,
            'observation_covariance': [[0.1, 0], [0, 0.1]]
        }))

    def test_make_kalman_filter(self):
        t = list(self.simple.tracks)[0]
        kf = t.make_kalman_filter(self.mp)
        self.assertFalse(kf is None)
        self.assertFalse(kf.predictions is None)
        self.assertFalse(kf.corrections is None)

        preds = list(kf.predictions)
        self.assertTrue(len(preds) > 0)
        corrs = list(kf.corrections)
        self.assertTrue(len(corrs) > 0)

        corr_states = list(x.state for x in corrs)
        corr_covs = list(x.covariance for x in corrs)
        self.assertEqual(len(corr_states), len(corrs))
        self.assertEqual(len(corr_covs), len(corrs))
        for s, c in zip(corr_states, corr_covs):
            self.assertEqual(len(s), 4)
            self.assertEqual(len(c), 4)
            self.assertEqual(len(c[0]), 4)
            self.assertEqual(len(c[1]), 4)
            self.assertEqual(len(c[2]), 4)
            self.assertEqual(len(c[3]), 4)

        try:
            import numpy as np
            corr_states = list(np.array(x.state) for x in corrs)
            corr_covs = list(np.array(x.covariance) for x in corrs)
        except ImportError:
            # it's ok not to have numpy
            pass

if __name__ == '__main__':
    unittest.main(argv=sys.argv[:1])
