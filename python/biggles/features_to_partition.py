"""Convert features.dat-style files from quincy projects to Biggles partitions."""

from __future__ import print_function

import argparse
import hashlib
import json
import logging
import struct
import sys

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
log = logging.getLogger()

def setup_parser(parser):
    parser.add_argument('input', metavar='FILE', type=str,
            help='Read features from FILE. If FILE is "-" read from standard input.')
    parser.add_argument('output', metavar='FILE', type=str,
            help='Write JSON formatted output to FILE. If FILE is "-" write to standard output.')

class TextFeatureParseError(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return self.msg

def read_text_features(input_file):
    lines = [l.strip() for l in input_file]
    if len(lines) == 0:
        raise TextFeatureParseError('Features file is empty')

    has_bg = False
    if lines[0].startswith('Flags:'):
        flags = lines.pop(0).split(':')[1:]
        has_bg = 'HasBg' in flags

    if len(lines) == 0:
        raise TextFeatureParseError('Features file has not frame count')
    n_frames = int(lines.pop(0))

    features = []
    for frame_idx in xrange(n_frames):
        if len(lines) == 0:
            raise TextFeatureParseError('Features file ends prematurely when reading frame %s' % (frame_idx,))

        frame_id, n_features = [int(t) for t in lines.pop(0).split(' ')]
        if len(lines) < n_features:
            raise TextFeatureParseError('Not enough lines in file when reading features for frame %s' % (frame_idx,))

        for feature_id in xrange(n_features):
            x, y, t = (0.0, 0.0, 0)

            fields = lines.pop(0).split(' ')

            has_pos = int(fields.pop(0))
            if has_pos != 0:
                x = float(fields.pop(0))
                y = float(fields.pop(0))

            # FIXME: interpolated, intensities, etc.
            features.append((x,y,frame_id))

    return features

def main(args):
    input_file = sys.stdin
    if args.input != '-':
        input_file = open(args.input, 'r')

    output_file = sys.stdout
    if args.output != '-':
        output_file = open(args.output, 'w')

    features = read_text_features(input_file)
    json.dump({ 'clutter': features, 'tracks': [] }, output_file)

    return 0

if __name__ == '__main__':
    # Parse the command line options
    parser = argparse.ArgumentParser(description=globals()['__doc__'])
    setup_parser(parser)
    args = parser.parse_args()

    # Run the main  program
    sys.exit(main(args))
