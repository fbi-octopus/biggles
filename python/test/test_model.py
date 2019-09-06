import json
import StringIO
import os
import unittest
import sys

from biggles import ModelParameters, Partition, ObservationCollection, TrackCollection

test_data_path = sys.argv[1]

class MPTestCase(unittest.TestCase):
    def setUp(self):
        self.mp = ModelParameters.from_json(json.dumps({
            'births_per_frame': 0.0,
            'clutter_per_frame': 0.0,
            'survival_probability': 0.0,
            'observation_probability': 0.0,
            'observation_covariance': [[0.1, 0], [0, 0.1]]
    }))

    def test_default_values(self):
        self.assertEqual(self.mp.mean_new_tracks_per_frame, 0.0)
        self.assertEqual(self.mp.mean_false_observations_per_frame, 0.0)
        self.assertEqual(self.mp.frame_to_frame_survival_probability, 0.0)
        self.assertEqual(self.mp.generate_observation_probability, 0.0)

    def test_non_default_values(self):
        self.mp.mean_new_tracks_per_frame = 1.0
        self.mp.mean_false_observations_per_frame = 2.0
        self.mp.frame_to_frame_survival_probability = 3.0
        self.mp.generate_observation_probability = 4.0

        self.assertAlmostEqual(self.mp.mean_new_tracks_per_frame, 1.0)
        self.assertAlmostEqual(self.mp.mean_false_observations_per_frame, 2.0)
        self.assertAlmostEqual(self.mp.frame_to_frame_survival_probability, 3.0)
        self.assertAlmostEqual(self.mp.generate_observation_probability, 4.0)

    def test_json(self):
        d = json.loads(self.mp.to_json())
        self.assertTrue('observation_probability' in d)

        mp_json = json.dumps({
                'births_per_frame': 1.0,
                'clutter_per_frame': 1.0,
                'survival_probability': 0.9,
                'observation_probability': 0.9,
                'observation_covariance': [[0.1, 0], [0, 0.1]]
        })
        mp = ModelParameters.from_json(mp_json)
        mp_d = json.loads(mp.to_json())

        self.assertTrue('births_per_frame' in mp_d)
        self.assertEqual(mp_d['births_per_frame'], 1.0)

class PartitionTestCase(unittest.TestCase):
    def setUp(self):
        self.p = Partition()

        p_json = open(os.path.join(test_data_path, 'simple-ground-truth.json')).read()
        self.simple = Partition.from_json(p_json)

    def test_default_values(self):
        self.assertEqual(self.p.first_time_stamp, 0)
        self.assertEqual(self.p.last_time_stamp, 0)
        self.assertEqual(self.p.duration, 0)
        self.assertTrue(isinstance(self.p.clutter, ObservationCollection))
        self.assertTrue(isinstance(self.p.tracks, TrackCollection))

    def test_json(self):
        d = json.loads(self.p.to_json())
        self.assertTrue('clutter' in d)
        self.assertTrue('tracks' in d)

        p_json = open(os.path.join(test_data_path, 'simple-ground-truth.json')).read()
        p_good = json.loads(p_json)
        p = Partition.from_json(p_json)
        p_d = json.loads(p.to_json())

        self.assertTrue('clutter' in p_d)
        self.assertEqual(len(p_d['clutter']), len(p_good['clutter']))
        self.assertTrue('tracks' in p_d)
        self.assertEqual(len(p_d['tracks']), len(p_good['tracks']))

    def test_clutter(self):
        clutter = self.simple.clutter
        self.assertNotEqual(len(clutter), 0)
        self.assertTrue(clutter.last_time_stamp > clutter.first_time_stamp)

        last_t = clutter.first_time_stamp - 1
        for idx, p in enumerate(clutter):
            self.assertTrue(p.x >= 0.0)
            self.assertTrue(p.y >= 0.0)
            self.assertTrue(p.t >= clutter.first_time_stamp)
            self.assertTrue(p.t < clutter.last_time_stamp)
            self.assertTrue(p.t >= last_t)

            last_t = p.t
        self.assertEqual(idx, len(clutter)-1)

    def test_Tracks(self):
        tracks = self.simple.tracks
        self.assertNotEqual(len(tracks), 0)

        for idx, t in enumerate(tracks):
            self.assertTrue(t.duration > 0)
            self.assertTrue(len(t.observations) > 0)

            self.assertTrue(t.first_time_stamp >= 0)
            self.assertTrue(t.last_time_stamp >= t.first_time_stamp)

            last_t = t.observations.first_time_stamp - 1
            for obs_idx, p in enumerate(t.observations):
                self.assertTrue(p.x >= 0.0)
                self.assertTrue(p.y >= 0.0)
                self.assertTrue(p.t >= t.observations.first_time_stamp)
                self.assertTrue(p.t < t.observations.last_time_stamp)
                self.assertTrue(p.t >= last_t)

                last_t = p.t
            self.assertEqual(obs_idx, len(t.observations)-1)
        self.assertEqual(idx, len(tracks)-1)

if __name__ == '__main__':
    unittest.main(argv=sys.argv[:1])
