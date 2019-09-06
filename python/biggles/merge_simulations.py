# -*- coding: utf-8 -*-
"""
Merging the data from the second input file into the data of the first input file.

Created on Wed Mar 20 10:25:12 2019

@author: vcn81216
"""

import argparse
import sys
#import os
import json
import numpy as np

def setup_parser(parser):
    parser.add_argument("files", nargs=2, help='biggles simulation files', metavar='file')
    parser.add_argument("-o", "--output", help="combine the input files into a single file", type=str)
    parser.add_argument("-v", "--verbose", action = 'store_true', default = False, help = "verbose")


class Merger(object):
    def __init__(self, args):
        self._args = args

    def get_time_stamps(self, tracks):
        time_stamps = set()
        for track in tracks:
            for obs in track['observations']:
                time_stamps.add(int(obs[2]))
        return time_stamps

    def load_input(self, fname):
        with open(fname) as fh:
            input_partition = json.load(fh)
        return input_partition

    def dump_partition(self, output_partition):
        if self._args.output is None:
            print(json.dumps(output_partition))
        else:
            with open(self._args.output, 'w') as fh:
                json.dump(output_partition, fh)

    def offset_tracks(self, tracks, offset):
        incr = lambda x : x + offset
        for track in tracks:
            track['observations'] = map(incr, track['observations'])
        return tracks

    def main(self):
        fn1, fn2 = self._args.files
        partition = self.load_input(fn1)
        addition = self.load_input(fn2)
        num_obs = len(partition['observations'])
        partition['observations'].extend(addition['observations'])
        partition['ground_truth'].extend(addition['ground_truth'])
        partition['unobserved'].extend(addition['unobserved'])
        partition['initial_partition']['tracks'].extend(
            self.offset_tracks(addition['initial_partition']['tracks'], num_obs))
        partition['initial_partition']['clutter'].extend(
            map(lambda x : x + num_obs, addition['initial_partition']['clutter']))
        partition['tracks'].extend( self.offset_tracks(addition['tracks'], num_obs) )
        partition['clutter'].extend( map(lambda x : x + num_obs, addition['clutter']) )
        partition['metadata']['merged_data'] = {
            'metadata' : addition['metadata'],
            'fileheader' : addition['fileheader'],
            'simulation' : addition['simulation'],
        }

        self.dump_partition(partition)

def main(args):
    run = Merger(args)
    run.main()

if __name__ == '__main__':
    # Parse the command line options
    parser = argparse.ArgumentParser(description=globals()['__doc__'],
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    setup_parser(parser)
    args = parser.parse_args()
    # Run the main  program
    sys.exit(main(args))
