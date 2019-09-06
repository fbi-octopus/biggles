import sys
import unittest

from biggles import TrackingState, Engine, ModelParameters, Partition, MoveType
from biggles import Observation
from biggles import ObservationCollection
from biggles import TrackCollection

class TSTestCase(unittest.TestCase):
    def setUp(self):
        self.ts = TrackingState()

    def test_attributes(self):
        self.assertTrue(hasattr(self.ts, 'sample_count'))
        self.assertTrue(hasattr(self.ts, 'current_log_pdf'))
        self.assertTrue(hasattr(self.ts, 'current_model_parameters'))
        #self.assertTrue(hasattr(self.ts, 'current_partition'))
        self.assertTrue(hasattr(self.ts, 'best_log_pdf'))
        self.assertTrue(hasattr(self.ts, 'best_model_parameters'))
        #self.assertTrue(hasattr(self.ts, 'best_partition'))
        self.assertTrue(hasattr(self.ts, 'acceptance_rate'))
        self.assertTrue(hasattr(self.ts, 'move_histogram'))
        self.assertTrue(hasattr(self.ts, 'current_move_type'))

    def test_default_values(self):
        self.assertEqual(self.ts.sample_count, 0)
        self.assertEqual(self.ts.current_log_pdf, 0.0)
        self.assertEqual(self.ts.best_log_pdf, 0.0)
        self.assertEqual(self.ts.acceptance_rate, 0.0)
        self.assertEqual(self.ts.current_move_type, MoveType.NONE)
        self.assertEqual(len(self.ts.move_histogram), 0)
        for v in self.ts.move_histogram:
            self.assertEqual(v, 0)

    def test_invalid_values(self):
        def bad():
            self.ts.sample_count = 'hello'
            self.ts.current_log_pdf = 'hello'
            self.ts.best_log_pdf = 'hello'
            self.ts.acceptance_rate = 'hello'

        self.assertRaises(Exception, bad)

    def test_non_default_values(self):
        self.ts.sample_count = 5
        self.ts.current_log_pdf = 1.3
        self.ts.best_log_pdf = 8.5
        self.ts.acceptance_rate = 1.9

        self.assertAlmostEqual(self.ts.sample_count, 5)
        self.assertAlmostEqual(self.ts.current_log_pdf, 1.3)
        self.assertAlmostEqual(self.ts.best_log_pdf, 8.5)
        self.assertAlmostEqual(self.ts.acceptance_rate, 1.9)

class EngineTestCase(unittest.TestCase):
    def setUp(self):
        self.e = Engine()
        oc = []
        oc.append(Observation.new(1,0,1))
        oc.append(Observation.new(0,1,2))
        oc.append(Observation.new(0,0,3))
        tc = TrackCollection()
        mp = ModelParameters()
        mp.frame_to_frame_survival_probability = 0.5
        mp.generate_observation_probability = 0.5
        mp.mean_false_observations_per_frame = 0.01
        mp.mean_new_tracks_per_frame = 0.01
        mp.initQ([0.1, 0.01]*2)
        obs_cov = [[1.0, 0.0], [0.0, 1.0]]
        mp.observation_error_covariance = obs_cov
        self.e.set_initial_conditions(mp, Partition.new(tc, oc))

    #def test_initial_conditions(self):
    #    oc = ObservationCollection()
    #    oc.insert(Observation.new(0,0,1))
    #    oc.insert(Observation.new(0,0,2))
    #    tc = TrackCollection()
    #    self.e.set_initial_conditions(ModelParameters(), Partition.new(tc, oc))


    def test_start_stop(self):
        import time
        self.assertEqual(self.e.tracking_state.sample_count, 0)
        self.e.start()
        while self.e.tracking_state.sample_count < 10:
            time.sleep(0.05)
        self.e.stop()

if __name__ == '__main__':
    unittest.main(argv=sys.argv[:1])
