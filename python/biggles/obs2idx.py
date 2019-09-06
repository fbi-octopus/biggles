# -*- coding: utf-8 -*-
"""
Convert observation partition into indexed partition

Created on Wed Sep 21 09:31:47 2016
@author: vcn81216
"""

from __future__ import print_function

import argparse
import sys
from indexing import has_flag, Flag, indexed_partition
import json

def setup_parser(parser):
    '''
    command line argument parser
    '''
    parser.add_argument("files", nargs='*', help='json files to convert', metavar='file')
    
def convert2idx(fname):
    fidx = fname.replace('.input.json', '.idx.input.json')
    if fidx == fname:
        print("Couldn't change file name. Do nothing.")
        return
    with open(fname) as fh:
        data = json.load(fh)
    if has_flag(data, Flag.INDEXED):
        return
    with open(fidx, 'w') as fh:
        json.dump(indexed_partition(data), fh)
    print('Wrote "%s"' % fidx)
    
def main(args):
    for fname in args.files:
        convert2idx(fname)

if __name__ == '__main__':
    # Parse the command line options
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description=globals()['__doc__'])
    setup_parser(parser)
    args = parser.parse_args()
    # Run the main  program
    sys.exit(main(args))
