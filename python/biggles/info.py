#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Prints info about the input files.

Created on Tue Mar 29 12:07:43 2016

@author: vcn81216
"""

import sys
import argparse
import json
import numpy as np
from collections import defaultdict

from biggles._chunking import DataChunk
from biggles.ripley_k import ripley_l

def setup_parser(parser):
    """Create an argparse parser for the command line options."""

    parser.add_argument('input', metavar='INPUT', type=str, nargs='*',
        help='read tracking result from file named INPUT')
    parser.add_argument('-x', nargs = 2, type=int, help='x-range')
    parser.add_argument('-y', nargs = 2, type=int, help='y-range')
    parser.add_argument('-t', nargs = 2, type=int, help='t-range')
    parser.add_argument('--count', type = int, default = 0, 
                        help = 'target number of observation per chunk')
    parser.add_argument('--apron', type = int, default = 4, 
                        help = 'size of the chunk overlap')
        
def get_range(amin, amax, lo, hi):
    assert(lo < hi)
    return (max(amin, lo), min(amax, hi))
    
def observe_parameters(data):
    clutter = np.array(data['clutter'])
    xmin, ymin, tmin = np.min(clutter, axis=0)
    xmax, ymax, tmax = np.max(clutter, axis=0)
    #clutter_count = dict()
    #for t in range(tmin, tmax+1):
    #    clutter_count[t] = len([x for x in clutter if x[2] == t])
    tracks = data['tracks']
    num_tracks  = len(tracks)
    volume = float((xmax - xmin) * (ymax - ymin))
    factor =  10000.0/volume
    if num_tracks == 0:
        return (0.0, factor*float(len(clutter))/float(tmax-tmin+1), 0, 0)
    total_length = 0
    obs_in_tracks = 0
    for t in tracks:
        obs_in_tracks += len(t['observations'])
        total_length += t['time_span'][1] - t['time_span'][0] + 1
        tmin = min(tmin, t['time_span'][0])
        tmax = max(tmax, t['time_span'][1])
        xmin0, ymin0, _ = np.min(t['observations'], axis=0)
        xmax0, ymax0, _ = np.max(t['observations'], axis=0)
        xmin, xmax = min(xmin0, xmin), max(xmax0, xmax)
        ymin, ymax = min(ymin0, ymin), max(ymax0, ymax)
    duration = float(tmax-tmin+1)
    volume = float((xmax - xmin) * (ymax - ymin))
    factor =  10000.0/volume
    #print("volume = %.1f" % volume)
    #print("duration = %d" % duration)
    #print("num_tracks = %d" % num_tracks)
    #print("clutter size = %d" % len(clutter))
    return (
        float(num_tracks)/duration * factor,
        float(len(clutter))/duration * factor,
        float(obs_in_tracks)/float(total_length),
        1.0 - float(num_tracks)/float(total_length))
    
def show_info(data, args):
    '''print some info about the biggles input data'''
    clutter = np.array(data['clutter'])
    xmin, ymin, tmin = np.floor(np.min(clutter, axis=0)).astype(int)
    xmax, ymax, tmax = np.ceil(np.max(clutter, axis=0)).astype(int)
    print('clutter size = {0}'.format(len(clutter)))
    print(' [x, y, t] = [[{0}, {1}], [{2}, {3}], [{4}, {5}]]'.format(
        xmin, xmax, ymin, ymax, tmin, tmax
        ))
    if False:
        ripley = ripley_l(clutter, 0.9, np.arange(.25, 6.1, 0.25))
        if len(ripley):
            print("Ripley's L function:")
            for key in sorted(ripley.keys()):
                print("* {0:.2f} => {1:.2f}".format(key, ripley[key]-key))
        print("b = {0:.4f}, c = {1:.4f}, o = {2:.4f}, s = {3:.4f}".format(*observe_parameters(data)))
    limited = False
    if not args.x is None:
        xmin, xmax = get_range(xmin, xmax, *args.x)
        limited = True
    if not args.y is None:
        ymin, ymax = get_range(ymin, ymax, *args.y)
        limited = True
    if not args.t is None:
        tmin, tmax = get_range(tmin, tmax, *args.t)
        limited = True
    if limited:
        rmin = np.array([xmin, ymin, tmin])
        rmax = np.array([xmax, ymax, tmax])
        selected = [ obs for obs in clutter if (rmin <= obs).all() and (obs <= rmax).all() ]
        print('selected size = {0}'.format(len(selected)))
        print(' [x, y, t] = [[{0}, {1}], [{2}, {3}], [{4}, {5}]]'.format(
            xmin, xmax, ymin, ymax, tmin, tmax
            ))

def show_chunk(data, args):
    '''print the chunk size for a desired number of observation per chunk'''
    chunklist = list()
    chunk0 = DataChunk.from_partition(data)
    chunk0.set_apron(args.apron)
    chunk0.set_size(args.count)
    
    print(chunk0)
    chunk0.get_chunks(chunklist)
    #get_chunks(chunk0, target, 0, chunklist)
    #header = '{0:8s} {1:4s} {2:4s} {3:16s} {4:16s}'.format('chunk id', 'row', 'col', 'x-range', 'y-range')
    #print(header)
    areas = defaultdict(int)
    for n, chunk in enumerate(chunklist):
        areas[int(chunk.area())] += 1
        print("{0:3d} {1:s}".format(n+1, chunk))
    print('{0:>6s}: {1:s}'.format('area', 'freq'))
    for area in sorted(areas.keys()):
        print('{0:6d}: {1:3d}'.format(area, areas[area]))
        
def main(args):
    '''the main function call'''
    for fname in args.input:
        with open(fname) as fh:
            data = json.load(fh)
        if args.count > 0:
            show_chunk(data, args)
        else:
            show_info(data, args)

if __name__ == '__main__':
    # Parse the command line options
    parser = argparse.ArgumentParser(description=globals()['__doc__'],
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    setup_parser(parser)
    args = parser.parse_args()

    # Run the main  program
    sys.exit(main(args))
