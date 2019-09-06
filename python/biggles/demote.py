"""Demote observations from tracks in a Biggles partition to the clutter."""

from __future__ import print_function

import argparse
import json
import logging
import sys

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
log = logging.getLogger()

def setup_parser(parser):
    parser.add_argument('input', metavar='INPUT', type=str, default='-', nargs='?',
            help='Write JSON formatted output to INPUT. If INPUT is "-" read from standard output.')
    parser.add_argument('output', metavar='OUTPUT', type=str, default='-', nargs='?',
            help='Write JSON formatted output to OUTPUT. If OUTPUT is "-" write to standard output.')

def main(args):
    input_file = sys.stdin
    if args.input is None or args.input != '-':
        input_file = open(args.input, 'r')

    output_file = sys.stdout
    if args.output is None or args.output != '-':
        output_file = open(args.output, 'w')

    partition = json.load(input_file)
    tracks = partition['tracks']
    clutter = partition['clutter']
    partition['initial_partition'] = { 'tracks' : list(tracks), 'clutter' : list(clutter) }

    for track in tracks:
        for obs in track['observations']:
            clutter.append(obs)
    partition['tracks'] = []

    json.dump(partition, output_file)

    return 0

if __name__ == '__main__':
    # Parse the command line options
    parser = argparse.ArgumentParser(description=globals()['__doc__'])
    setup_parser(parser)
    args = parser.parse_args()

    # Run the main  program
    sys.exit(main(args))

# vim:tw=120
