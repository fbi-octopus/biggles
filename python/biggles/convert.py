"""Various CLI tools to convert file formats to/from Biggles partition formats."""

import argparse
from collections import deque
import json
import sys
import os
import numpy as np
from _filedescriptor import filedescriptor

class Feature(object):
    ''' a feature that has a position and intensities in one or more channels
    '''
    def __init__(self):
        self._x = None
        self._y = None
        self._evidence = None
        self._interp = None
        self._sigma = None
        self._frame = None
        self._channels = {}
        
    def __str__(self):
        return "Feature (x={0}, y={1}, t={2})".format(self.x(), self.y(), self.frame())
        
    def from_parser(self, frame, feat, records):
        self._frame = frame
        self._x = feat['x']
        self._y = feat['y']
        self._evidence = feat['evidence']
        self._interp = feat['interp']
        self._sigma = None
        if 'sgima' in feat:
            self._sigma = feat['sigma']
        ''' "counts", "sigcounts", "bg" '''
        self._channels = records
        return self
        
    def close_to(self, other, epsilon):
        return (np.hypot(self.x()-other.x(), self.y()-other.y()) < epsilon
            and other.frame == self.frame)
        
    def __iter__(self):
        return self._channels.__iter__()
        
    def __getitem__(self, key):
        return self._channels[key]
        
    def x(self):
        return self._x
    def y(self):
        return self._y
    def loc(self):
        return (self.x(), self.y())
    def channels(self):
        return self._channels.keys()
    def frame(self):
        return self._frame
    def is_interpolated(self):
        return self._interp
    def is_detected(self):
        ''' was feature detected '''
        return not (self._interp is None or self._interp)
        
    def obs(self):
        ''' Returns the feature as a Biggles style observation
            Return (x, y, t), where 't' is the frame id
        '''
        return (self.x(), self.y(), self.frame())
        
    def __cmp__(self, other):
        if self.frame() > other.frame():
            return 1
        if self.frame() < other.frame():
            return -1
        if self.x() > other.x():
            return 1
        if self.x() < other.x():
            return -1
        return (self.y() > other.y()) - (self.y() < other.y())
        
#    def __lt__(self, other):
#        return self.x() < other.x() or self.x() == other.x() and self.y() < other.y()
#    def __eq__(self, other):
#        return self.x() == other.x() and self.y() == other.y()
#    def __ne__(self, other):
#        return self.x() != other.x() or self.y() != other.y()
#    def __gt__(self, other):
#        return self.x() > other.x() or self.x() == other.x() and self.y() > other.y()
#    def __le__(self, other):
#        return self.x() < other.x() or self.x() == other.x() and self.y() <= other.y()
#    def __ge__(self, other):
#        return self.x() > other.x() or self.x() == other.x() and self.y() >= other.y()
        
class Track(object):
    def __init__(self, id=None, size=None):
        self.set_id(id)
        self._features = []
        self._size = size
        
    def __len__(self):
        return len(self._features)
    def __getitem__(self, key):
        return self._features[key]
    def __iter__(self):
        return self._features.__iter__()
    
    def set_id(self, id):
        assert(id is None or isinstance(id, int) and id > 0)
        self._id = id
    
    def id(self):
        return self._id
        
    def append(self, feature):
        ''' appends a feature to the list of features
        '''
        assert(isinstance(feature, Feature))
        self.get().append(feature)
        
    def frames(self):
        ''' returns the *sorted* list of frames occuring in the track
        '''
        return sorted(set([ feat.frame() for feat in self.get()]))
        
    def frames_spanned(self):
        ''' Returns the number of frames that are covered by this track
            This can be larger than the lenght of the track (which is the number of features).
            In some frames the track might not have a feature
        '''
        frames = self.frames()
        return frames[-1] - frames[0] + 1
        

    def bounding_box(self):
        ''' returns a box ((x_min, y_min), (x_max, y_max)) that conatains all features
        '''
        xs = [ feat.x() for feat in self.get() ]
        ys = [ feat.y() for feat in self.get() ]
        return ((min(xs), min(ys)), (max(xs), max(ys)))
        
    def get(self):
        ''' returns the list of features as "Feature" objects
        '''
        return self._features
        
    def as_obs(self):
        ''' Return the list of features as Biggles style observation (x, y, t)
        '''
        return [ feat.obs() for feat in self.get() ]

def _parse_feature(fields, has_bg):
    feature = {}
    records = {}

    has_pos = int(fields.pop(0))
    if has_pos != 0:
        feature['x'] = float(fields.pop(0))
        feature['y'] = float(fields.pop(0))

    has_evidence = int(fields.pop(0))
    if has_evidence != 0:
        feature['evidence'] = float(fields.pop(0))

    #if int(fields.pop(0)) != 0:
    #    feature['interp'] = True
    feature['interp'] = int(fields.pop(0)) != 0

    n_channels = int(fields.pop(0))
    for ch_idx in range(n_channels):
        record = {}
        ch_id = int(fields.pop(0))
        record['counts'] = float(fields.pop(0))
        record['sigcounts'] = float(fields.pop(0))
        if has_bg:
            record['bg'] = float(fields.pop(0))
        records[ch_id] = record

    # Profile sigma
    if len(fields) > 0:
        hassigma = int(fields.pop(0)) != 0
        if hassigma:
            feature['sigma'] = float(fields.pop(0))

    return feature, records

def _read_text_features(input_file, channel_q=0):
    """Read text features from a MSMM features.dat file. The channel_q argument is the channel id to extract."""

    lines = deque(l.strip() for l in input_file)
    if len(lines) == 0:
        raise Exception('Features file is empty')

    has_bg = False
    if lines[0].startswith('Flags:'):
        flags = lines.popleft().split(':')[1:]
        has_bg = 'HasBg' in flags

    if len(lines) == 0:
        raise Exception('Features file has not frame count')
    n_frames = int(lines.popleft())

    features = []
    frame_ids = []
    for frame_idx in xrange(n_frames):
        if len(lines) == 0:
            raise Exception('Features file ends prematurely when reading frame %s' % (frame_idx,))

        frame_id, n_features = [int(t) for t in lines.popleft().split(' ')]
        frame_ids.append(frame_id)
        if len(lines) < n_features:
            raise Exception('Not enough lines in file when reading features for frame %s' % (frame_idx,))

        # basically a transliteration of the logic in def_Feature.hpp:Feature::Feature
        for idx in xrange(n_features):
            fields = lines.popleft().split(' ')
            feature, new_records = _parse_feature(fields, has_bg)

            if channel_q in new_records:
                r = new_records[channel_q]
                r['frame'] = frame_id
                r.update(feature)
                features.append(r)
            elif channel_q == 'all':
                features.append(Feature().from_parser(frame_id, feature, new_records))

    return (features, frame_ids)

def _parse_track_frame(line):
    ''' Converts a MSMM ".dat" style feature into a "Feature" object
    '''
    data = line.split(' ')
    frame_id = int(data.pop(0))
    return Feature().from_parser(frame_id, *_parse_feature(data, True))
    
def _load_simple_tracks(fname):
    ''' Loads "tracks_Simple.data"
        Returns = [ (track_id, [ feature, ...]), ...]
        If the data contains multi-channel features, then each feature occures in the track
        as frequently as the number of channels in which it appears
    '''
    with open(fname) as fh:
        lines = [l.strip() for l in fh]
    flags = lines.pop(0)
    assert(flags == "Flags:HasBg")
    num_tracks = int(lines.pop(0))
    tracks = []
    while lines:
        track_id, track_size = map(int, lines.pop(0).split(' '))
        track = Track(track_id, track_size)
        for frame in range(track_size):
            track.append(_parse_track_frame(lines.pop(0)))
        tracks.append(track)
    assert(len(tracks) == num_tracks)
    return tracks
        
def _time_chunk_features(features, frame_ids, n_chunks):
    '''
    returns a list of dictionaries that contain features, first frame in chunk, last frame in chunk
    '''
    n_frames = len(frame_ids)
    if n_frames < n_chunks:
        sys.stderr.write('Input file has not enough frames for {0} time chunks\n'.format(n_chunks))
        return None
    base_len = n_frames / n_chunks
    num_ext  = n_frames % n_chunks
    ext = [base_len]*n_chunks
    ext[0:num_ext] = [base_len+1]*num_ext
    hi = 0
    res = []
    for idx in xrange(n_chunks):
        this = {}
        lo = hi
        hi = lo + ext[idx]
        frame_range = frame_ids[lo:hi]
        this["features"] = [item for item in features if item["frame"] in frame_range]
        this["lo_frame"] = frame_range[0]
        this["hi_frame"] = frame_range[-1]
        res.append(this)
    return res


def _msmm_arg_parser(args=None):
    parser = argparse.ArgumentParser(description='Convert MSMM features.dat files to Biggles partitions.')
    parser.add_argument('features', metavar='FILE', type=str,
            help='MSMM features.dat file to convert or - for standard input')
    parser.add_argument('channel', metavar='CHANNEL', type=int,
            help='Channel id to select')
    parser.add_argument('--output', '-o', metavar='FILE', type=str,
            help='Write output to FILE instead of standard output')
    parser.add_argument('--time-chunks','-t', type=int, default = 1,
            help='Number of temporal chunks that will be created')

    return parser.parse_args(args)

def msmm_to_biggles(args=None):
    """Convert MSMM Features files to Biggles partitions. If args is None, use sys.argv, otherwise args is a sequence of
    strings to interpret as the command line arguments."""

    # Parse command line arguments
    args = _msmm_arg_parser(args)
    project_id = "stdin"

    input_file = sys.stdin
    if 'features' in args:
        input_file = open(args.features)
        project_id = os.path.basename(os.path.dirname(os.path.abspath(args.features)))


    # Reat in the features
    (features, frame_ids) = _read_text_features(input_file, args.channel)
    if len(features) == 0:
        sys.stderr.write('Input file has no features for channel {0}\n'.format(args.channel))
        sys.exit(1)
    if args.time_chunks > 1:
        chunked_features = _time_chunk_features(features, frame_ids, args.time_chunks)
        for entry in chunked_features:
            # determine filename and open file
            fname = "features.json"
            if 'output' in args and args.output:
                fname = args.output
            (base, ext) = os.path.splitext(fname)
            output_file = open(base + "_frames_" + str(entry["lo_frame"]) + "-" + str(entry["hi_frame"]) + ext, 'w')
            # Write output
            json.dump({
                'fileheader' : filedescriptor(),
                'metadata': {"project": project_id, "channel": args.channel},
                'tracks': [],
                'clutter': list((f['x'], f['y'], f['frame']) for f in entry["features"]),
            }, output_file)
    else:

        # Open the output file for writing
        output_file = sys.stdout
        if 'output' in args:
            output_file = open(args.output, 'w')

        # Write output
        json.dump({
            'fileheader' : filedescriptor(),
            'metadata': {"project": project_id, "channel": args.channel},
            'tracks': [],
            'clutter': list((f['x'], f['y'], f['frame']) for f in features),
        }, output_file)

    sys.exit(0)
