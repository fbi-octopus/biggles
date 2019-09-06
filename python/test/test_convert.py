import os
import sys
import unittest

from biggles.convert import msmm_to_biggles

test_data_path = sys.argv[1]

class MSMMConvertTestCase(unittest.TestCase):
    def setUp(self):
        self.features_dat_path = os.path.join(test_data_path, 'msmm-features.dat')
        self.output_path = os.path.join(test_data_path, 'msmm-features.json')

    def test_help_option(self):
        # Check that help causes a successful exit
        self.assertRaises(SystemExit, lambda: msmm_to_biggles(('-h',)))

    def test_features_exist(self):
        self.assertTrue(os.path.isfile(self.features_dat_path))

    def test_convert_invalid_channel(self):
        # Test that an invalid channel does not convert successfully
        for ch in (0, 4):
            self.assertRaises(SystemExit, lambda: msmm_to_biggles((self.features_dat_path, str(ch))))

    def test_convert_channel_1(self):
        # Test that channel 1 converts successfully
        self.assertRaises(SystemExit, lambda: msmm_to_biggles((self.features_dat_path, '1', '-o', self.output_path)))

    def test_convert_all_channels(self):
        # Test that channels 1, 2 and 3 convert successfully
        for ch in xrange(1, 4):
            self.assertRaises(SystemExit, lambda: msmm_to_biggles((self.features_dat_path, str(ch), '-o', self.output_path)))

if __name__ == '__main__':
    unittest.main(argv=sys.argv[:1])
