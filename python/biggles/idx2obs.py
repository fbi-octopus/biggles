# -*- coding: utf-8 -*-
"""
Convert indexed partition into observation partition

Created on Mon Jul 24 16:10:14 2017
@author: vcn81216
"""

from __future__ import print_function

import argparse
import sys
from indexing import has_flag, Flag, observation_partition
import json
import os

def setup_parser(parser):
    '''
    command line argument parser
    '''
    parser.add_argument("files", nargs='*', help='json files to convert', metavar='file')

def convert2obs(fname):
    fobs = fname.replace('.idx.input.json', '.input.json')
    fobs = fobs.replace('.input.json', '.obs.input.json')
    if fobs == fname:
        print("Couldn't change file name. Do nothing.")
        return
    if os.path.exists(fobs):
        print('File "{}" exists. Skipping.'.format(fobs))
        return
    with open(fname) as fh:
        data = json.load(fh)
    if not has_flag(data, Flag.INDEXED):
        return
    with open(fobs, 'w') as fh:
        json.dump(observation_partition(data), fh)
    print('Wrote "%s"' % fobs)

def main(args):
    for fname in args.files:
        convert2obs(fname)

if __name__ == '__main__':
    # Parse the command line options
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description=globals()['__doc__'])
    setup_parser(parser)
    args = parser.parse_args()
    # Run the main  program
    sys.exit(main(args))
