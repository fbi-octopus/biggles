# -*- coding: utf-8 -*-
"""
Calculates the GED bewteen two partitons

Created on Tue Jul 25 14:23:12 2017 @author: vcn81216
"""

import argparse
import sys
import json
import numpy as np
from collections import defaultdict
from _track_utils import convert_initial_partition
from os.path import basename

def setup_parser(parser):
    parser.add_argument("files", nargs=2, help='nano-positioning input files to test', metavar='file')
    #parser.add_argument("-o", "--output", help="combine the input files into a single file", type=str)
    parser.add_argument("-r", "--root", action = 'store_true', default = False,
                        help = "Use root instead of simulation as reference.")
    parser.add_argument("-t", "--tracked", action = 'store_true', default = False,
                        help = "use the second file name as key to store the results.")
    parser.add_argument('-o', '--output',
                        help='file to store the results')

def load_partition(fname):
    ''' simply load a partition from file and return it '''
    with open(fname) as fh:
            partition = json.load(fh)
    return partition

def dist(o1, o2):
    return np.sqrt((o1[0]-o2[0])**2 + (o1[1]-o2[1])**2)

def log10bin(val):
    if val==0:
        return -20
    return max(int(np.log10(val)), -20)

class Comparer(object):
    ''' compares a partition with other partitions with respect to the GED '''
    def __init__(self, args):
        self._args = args
        self._reference = None
        self._pool = None
        self._num_ref_links = None
        self._output = None

    def set_reference_partition(self):
        fname, _ = self._args.files
        partition = load_partition(fname)
        have_sims, self._reference, sims_partition, self._pool = convert_initial_partition(partition)
        if not self._args.root and not have_sims:
            raise Exception('No simulation found in {0}'.format(fname))
        if not self._args.root:
            self._reference = sims_partition
        self._num_ref_links = sum([
            len(track.observations) - 1 for track in self._reference.tracks
        ])

    def load_partition(self):
        _, fname = self._args.files
        return load_partition(fname)

    def compare_to_root(self, partition):
        ''' does a comaprision to the partition in "root" '''
        _, partition, _, _ = convert_initial_partition(partition, self._pool)
        return self._reference.ged(partition)

    def _comp_to_part_tracks(self, part_tracks):
        ''' calc the GED between the reference and the partition defined by part_tracks '''
        ''' part_tracks does not contain clutter. Is it needed? '''
        partition = { 'tracks' : [], 'clutter' : range(len(self._pool)) }
        for track_stub in part_tracks:
            partition['tracks'].append({'observations' : track_stub[0], 'time_span' : track_stub[1]})
            for idx in track_stub[0]:
                partition['clutter'].remove(idx)
        _, partition, _, _ = convert_initial_partition(partition, self._pool)
        return self._reference.ged(partition)

    def compare_to_sample(self, data, chunk=0):
        ''' does a comaprision to the sampled partition in "chunk" '''
        assert(0 <= chunk < len(data))
        ged_list = []
        for part_tracks in data[chunk]['track_samples']:
            ged_list.append(self._comp_to_part_tracks(part_tracks))
        return ged_list

    def read_output(self):
        if not self._args.output is None:
            data = {'root' : {}, 'samples' : {} }
            try:
                with open(self._args.output) as fh:
                    data = json.load(fh)
            except IOError as ex:
                print('Notice: "{0}"'.format(ex))
            finally:
                self._output = data

    def write_output(self):
        if not self._args.output is None:
            try:
                with open(self._args.output, 'w') as fh:
                    json.dump(self._output, fh)
            except Exception as ex:
                print(ex)

    def dump_ged(self, what, data):
        if not self._args.output is None:
            assert(what in ['root', 'samples'])
            key, tracked = self._args.files
            if self._args.tracked:
                key = basename(tracked)
            self._output[what][key] = data

    def compare_observations(self):
        f1, f2 = self._args.files
        obslist1 = load_partition(f1)['observations']
        obslist2 = load_partition(f2)['observations']
        time_off_count = 0
        dists = defaultdict(int)
        for o1, o2 in zip(obslist1, obslist2):
            if o2[2] != o1[2]:
                time_off_count += 1
            else:
                dists[log10bin(dist(o1, o2))] += 1
        print("Observation pool analysis:")
        print("number of time errors = {0}".format(time_off_count))
        print("distances:")
        for key, val in dists.items():
            print(" {0:3d} => {1:5d}".format(key, val))

    def print_root(self, ged):
        print("Root:")
        print("GED={0:6.2f}%".format(100.0*float(ged)/self._num_ref_links))

    def print_samples(self, ged_list):
        percs = np.percentile(ged_list, [0.0, 16, 50, 84, 100.0])
        print("Samples")
        print("{0:6.2f}%".format(100.0*float(percs[0])/self._num_ref_links))
        print("{0:6.2f}%".format(100.0*float(percs[1])/self._num_ref_links))
        print("{0:6.2f}%".format(100.0*float(percs[2])/self._num_ref_links))
        print("{0:6.2f}%".format(100.0*float(percs[3])/self._num_ref_links))
        print("{0:6.2f}%".format(100.0*float(percs[4])/self._num_ref_links))

    def main(self):
        self.set_reference_partition()
        self.compare_observations()
        assert(not self._pool is None)
        partition = self.load_partition()
        self.read_output()
        if 'tracks' in partition and 'clutter' in partition:
            ged = self.compare_to_root(partition)
            self.print_root(ged)
            self.dump_ged('root', ged/self._num_ref_links)
        if 'samples' in partition:
            ged_list = self.compare_to_sample(partition['samples'])
            self.print_samples(ged_list)
            self.dump_ged('samples', (np.array(ged_list)/self._num_ref_links).tolist())
        self.write_output()

def main(args):
    run = Comparer(args)
    run.main()

if __name__ == '__main__':
    # Parse the command line options
    parser = argparse.ArgumentParser(description=globals()['__doc__'],
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    setup_parser(parser)
    args = parser.parse_args()
    # Run the main  program
    sys.exit(main(args))
