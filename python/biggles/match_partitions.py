# -*- coding: utf-8 -*-
"""
Matching partition files via kD-tree. Given two files, an indexed partition file
and a plain partition file, match the observations in the second to the observations in the first
and create an indexed file.

Created on Thu Jul 20 10:25:12 2017

@author: vcn81216
"""

import argparse
import sys
#import os
import json
from indexing import Flag, has_flag, indexed_partition, replace_flag
from _kdtree import KDTree
from _filedescriptor import filedescriptor
import numpy as np
from collections import defaultdict

def setup_parser(parser):
    parser.add_argument("files", nargs=2, help='nano-positioning input files to test', metavar='file')
    parser.add_argument("-o", "--output", help="combine the input files into a single file", type=str)
    parser.add_argument("-v", "--verbose", action = 'store_true', default = False, help = "verbose")

def _dist(p1, p2):
    return np.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)


class Matcher(object):
    def __init__(self, args):
        self._args = args
        self._trees = None
        self._observations = None
        self._indices_avail = None
        self._min_time_stamp = None
        self._num_obs_in_tracks = None
        self._verbose = self._args.verbose
        self._reference_partition = None
        self._input_partition = None

    def _setup_tree(self):
        ''' sets up the kD-trees and retrieves the min time stamp '''
        if not self._observations is None:
            timed_obs = defaultdict(list)
            for obs in self._observations:
                timed_obs[int(obs[2])].append(obs)
            self._trees = dict()
            for key, val in timed_obs.items():
                self._trees[key] = KDTree(val)
            self._min_time_stamp = min(timed_obs.keys())
        else:
            raise Exception("observations not set.")

    def make_all_obs_avail(self):
        self._indices_avail = range(len(self._observations))
        self._num_obs_in_tracks = 0

    def set_observations(self):
        self._observations = self._reference_partition['observations']
        self.make_all_obs_avail()

    def get_time_stamps(self, tracks):
        time_stamps = set()
        for track in tracks:
            for obs in track['observations']:
                time_stamps.add(int(obs[2]))
        return time_stamps

    def _track_obs_to_index(self, list_o_points):
        ''' for a list of points associate each point with an observation that is still available,
        i.e. which is still listed in the indices_avail
        '''
        indices = []
        for point in list_o_points:
            if np.isnan(point).any():
                continue
            point_copy = point[:]
            point_copy[2] += self._min_time_stamp
            obs = self._trees[int(point_copy[2])].query(point_copy)[0]
            if not (_dist(obs, point) < 0.001):
                print("minimal time stamp".format(self._min_time_stamp))
                print("Result of kD-tree query = {0}".format(obs))
                print("original point = {0}".format(point))
                assert(_dist(obs, point) < 0.001)
            obs_index = self._observations.index(obs)
            assert(obs_index in self._indices_avail)
            indices.append(obs_index)
            self._indices_avail.remove(obs_index)
        self._num_obs_in_tracks += len(indices)
        if self._verbose:
            print("#obs = {0}".format(len(indices)))
        return indices

    def time_span_from_observations(self, list_o_points):
        ''' getting the track time span from a list of points. It is defined as the pair of the
        time stamp of the earliest observation and the time stamp immediatly after the latest
        observation
        '''
        time_points = [ o[2] for o in list_o_points ]
        return [ int(min(time_points) + self._min_time_stamp ), int(max(time_points) + self._min_time_stamp + 1) ]

    def convert_tracks(self, tracks):

        return [{ 'observations' :  self._track_obs_to_index(track['observations']),
                  'time_span' : self.time_span_from_observations(track['observations'])
                } for track in tracks]

    def convert_partition(self):
        new_partition = {
            'observations' : self._observations,
            'tracks' : self.convert_tracks(self._input_partition['tracks']),
            'clutter' : self._indices_avail
        }
        ''' this doubling of the data is done to enable the 3D view in "biggles_vis_posterior.py"
        '''
        #new_partition['initial_partition'] = {
        #    'tracks' : new_partition['tracks'],
        #    'clutter' : new_partition['clutter']
        #}
        new_partition['initial_partition'] = self._reference_partition['initial_partition']
        replace_flag(new_partition, Flag.OBSERVATIONS, Flag.INDEXED)
        return new_partition

    def write_meta(self, output_partition):
        ''' writes the meta data to the output partition '''
        metadata = {}
        if 'simulation' in self._reference_partition:
            metadata['simulation'] = self._reference_partition['simulation']
        if 'fileheader' in self._reference_partition:
            metadata['sourceheader'] = self._reference_partition['fileheader']
        output_partition['metadata'] = metadata
        output_partition['fileheader'] = filedescriptor()

    def load_reference(self):
        reference_file, _ = self._args.files
        with open(reference_file) as fh:
            self._reference_partition = json.load(fh)
        if has_flag(self._reference_partition, Flag.INDEXED):
            pass
        else:
            self._reference_partition = indexed_partition(self._reference_partition)

    def load_input(self):
        _, input_file = self._args.files
        with open(input_file) as fh:
            self._input_partition = json.load(fh)

    def set_input_partition(self, partition):
        self._input_partition = partition

    def dump_partition(self, output_partition):
        if self._args.output is None:
            print(json.dumps(output_partition))
        else:
            with open(self._args.output, 'w') as fh:
                json.dump(output_partition, fh)

    def main(self):
        self.load_reference()
        self.load_input()

        self.set_observations()
        self._setup_tree()
        #print(self.get_time_stamps(input_partition['tracks']))

        output_partition = self.convert_partition()
        self.write_meta(output_partition)

        if self._verbose:
            print("total: {0}".format(self._num_obs_in_tracks))

        self.dump_partition(output_partition)

def main(args):
    run = Matcher(args)
    run.main()

if __name__ == '__main__':
    # Parse the command line options
    parser = argparse.ArgumentParser(description=globals()['__doc__'],
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    setup_parser(parser)
    args = parser.parse_args()
    # Run the main  program
    sys.exit(main(args))
