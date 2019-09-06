# -*- coding: utf-8 -*-
"""
Biggles sample recorder
Created on Thu Apr 14 15:08:32 2016

@author: vcn81216
"""

from collections import defaultdict
from _utilities import fifo
from _filedescriptor import filedescriptor
from numpy import floor, ceil, array
from _track_utils import extract_kalman_states, update_idx_pair_dicts
import json
import os
from datetime import datetime
from biggles import IndexedPartition
import hdf5util as h5u
import time
from itertools import izip

def extend_from_observations(observations):
    ''' return [xmin, xmax], [ymin, ymax] that contains all observations '''
    xx = [ o[0] for o in observations ]
    yy = [ o[1] for o in observations ]
    return [floor(min(xx)), ceil(max(xx))], [floor(min(yy)), ceil(max(yy))]

def samples_from_raw(raw_samples):
    '''converts raw samples into the output form'''
    samples = {
        "log_pdf" : list(raw_samples['log_pdf']),
        "sample_count" : raw_samples['sample_count'],
        'burnin_length' : raw_samples['burnin_length'],
        'acceptance_rate' : list(raw_samples['acceptance_rate']),
        'track_samples' : list(raw_samples['track_samples']),
        'chunk_id' : 0,
        'grid_row' : 0,
        'grid_column' : 0,
        'x_extent' : raw_samples['x_extent'],
        'y_extent' : raw_samples['y_extent'],
        'record_ticks' : raw_samples['record_ticks'],
    }
    if 'ged_to_sims' in raw_samples:
        samples['ged_to_sims'] = list(raw_samples['ged_to_sims'])
    samples['param_samples'] = {}
    for key, value in raw_samples['param_samples'].items():
        samples['param_samples'][key] = list(value)
    return samples

def get_modal_state(state):
    indexed_partition = IndexedPartition()
    indexed_partition.from_partition(state.best_partition)
    return {
        'log_pdf': state.best_log_pdf,
        'tracks': indexed_partition.tracks(),
        'clutter': indexed_partition.clutter(),
        'model_parameters': json.loads(state.best_model_parameters.to_json()),
    }

def partition_to_index_form(part):
    '''converts a partition from observation coordinates form to observation index form
    clutter is omitted (since it is 'observations' minus 'observations in tracks')
    uses the biggles.IndexedPartition class.
    A final track transformation is necessary (bummer!) that tranforms the dictionary form into
    the list form
    '''
    indexed_partition = IndexedPartition()
    indexed_partition.from_partition(part)
    index_form = []
    for track in indexed_partition.tracks():
        index_form.append([track['observations'], track['time_span']])
    return index_form

class SampleRecorder(object):
    ''' Records samples from biggles tracking and writes them.
    '''
    def __init__(self, initial_partition):
        self.initial_partition = initial_partition
        x_extent, y_extent = extend_from_observations(initial_partition['observations'])

        self.records = {
            "log_pdf" : list(),
            "sample_count" : 0,
            'burnin_length' : 0,
            'acceptance_rate' : list(),
            'param_samples' : defaultdict(list),
            'track_samples' : list(),
            'x_extent' : x_extent,
            'y_extent' : y_extent,
            'record_ticks' : 0,
        }
        self.biggles_version = 'N/A'
        self.last_state = None
        self.sims_partition = None

    def set_biggles_version(self, version):
        self.biggles_version = version

    def set_sims_partition(self, sims_partition):
        '''sets the realisation of the simulated partition'''
        self.sims_partition = sims_partition
        self.records['ged_to_sims'] = list()

    def record_state(self, state):
        '''records the sampling state'''
        if state.sample_count == 0:
            return;
        self.records['log_pdf'].append(state.current_log_pdf)
        if self.sims_partition is not None:
            self.records['ged_to_sims'].append(self.sims_partition.ged(state.current_partition))
        self.records['sample_count'] = state.sample_count
        self.records['acceptance_rate'].append(state.recent_acceptance_rate)
        parameters = json.loads(state.current_model_parameters.to_json())
        self.records['param_samples']['observation_covariance'].append(
            parameters['observation_covariance'])
        self.records['param_samples']['births_per_frame'].append(
            parameters['births_per_frame'])
        self.records['param_samples']['clutter_per_frame'].append(
            parameters['clutter_per_frame'])
        self.records['param_samples']['observation_probability'].append(
            parameters['observation_probability'])
        self.records['param_samples']['survival_probability'].append(
            parameters['survival_probability'])
        self.records['track_samples'].append(partition_to_index_form(state.current_partition))
        self.records['record_ticks'] += 1
        self.last_state = state

    def set_burnin_length(self, sample_count):
        ''' at the end of the burn-in phase the burnin length can be recorded by calling this. '''
        self.records['burnin_length'] = sample_count

    def write_hdf5_results(self, args, root):
        '''writes results to hdf5 handle'''
        metadata = {}
        if "metadata" in self.initial_partition:
            metadata.update(self.initial_partition["metadata"])
        if "fileheader" in self.initial_partition:
            metadata["sourceheader"] = self.initial_partition["fileheader"]
        if "simulation" in self.initial_partition:
            metadata["simulation"] = self.initial_partition["simulation"]
        if args.input is not None:
            metadata["inputfile"] = args.input
        descr = 'Biggles output. Either samples or averaged and not averaged link statistics'
        metadata["description"] = descr
        metadata['file version'] = 'biggles-out {0}'.format(self.biggles_version)
        root.attrs['metadata'] = json.dumps(metadata)


        if "ground_truth" in self.initial_partition:
            h5u.truth2hdf5(self.initial_partition['ground_truth'], root.create_group('ground_truth'))
        if "unobserved" in self.initial_partition:
            h5u.truth2hdf5(self.initial_partition['unobserved'], root.create_group('unobserved'))
        if "initial_partition" in self.initial_partition:
            h5u.write_hdf5_indexed_partition(self.initial_partition["initial_partition"], root)
        root.attrs["fileheader"] = json.dumps(filedescriptor())

        h5call = root.create_group("scalar_call_parameters")
        for key, val in vars(args).items():
            if isinstance(val, (int, float, bool, str)):
                h5call.attrs[key] = val
            if isinstance(val, type(None)):
                h5call.attrs[key] = json.dumps(val)
        root.create_dataset('observations', data = array(
            self.initial_partition['observations'], dtype = float))
        h5map = root.create_group('map_solution')
        h5map.attrs['log_pdf'] = self.last_state.best_log_pdf
        h5u.write_hdf5_parameters(self.last_state.best_model_parameters,
                                  h5map.create_group('model_parameters'))
        indexed_part = IndexedPartition()
        indexed_part.from_partition(self.last_state.best_partition)
        h5u.write_hdf5_indexed_partition(indexed_part, h5map)
        h5sam = root.create_group('samples')
        h5u.write_hdf5_samples(samples_from_raw(self.records), h5sam.create_group('0'))

    def write_results(self, args, output_file):
        '''writes all results to a json file'''
        output = {'metadata' : {}}
        if "metadata" in self.initial_partition:
            output["metadata"] = self.initial_partition["metadata"]
        if "fileheader" in self.initial_partition:
            output["metadata"]["sourceheader"] = self.initial_partition["fileheader"]
        if "simulation" in self.initial_partition:
            output["metadata"]["simulation"] = self.initial_partition["simulation"]
        if "ground_truth" in self.initial_partition:
            output["ground_truth"] = self.initial_partition['ground_truth']
        if "unobserved" in self.initial_partition:
            output["unobserved"] = self.initial_partition['unobserved']
        if "initial_partition" in self.initial_partition:
            output["initial_partition"] = self.initial_partition["initial_partition"]
        if args.input is not None:
            output["metadata"]["inputfile"] = args.input
        descr = 'Biggles output. Either samples or averaged and not averaged link statistics'
        output["metadata"]["description"] = descr
        output['metadata']['file version'] = 'biggles-out {0}'.format(self.biggles_version)
        output["fileheader"] = filedescriptor()
        paras = {}
        for key, val in vars(args).items():
            if isinstance(val, (str, int, float, bool, type(None), list, tuple)):
                paras[key] = val
        output["scalar_call_parameters"] = paras
        output['observations'] = self.initial_partition['observations']
        output['map_solution'] = get_modal_state(self.last_state)
        output['samples'] = [samples_from_raw(self.records)]
        json.dump(output, output_file)

class MpiSampleRecorder(object):
    '''Sample recorder for mpi_track'''
    def __init__(self, start_time):
        self.pair_dict = {}
        self.avg_pair_dict = dict()
        self.param_dict = {}
        self.param_dict['mean_track_length'] = []
        self.param_dict['track_count'] = []
        self.move_names = dict(zip(range(11), [
            'move_none', 'move_birth', 'move_death', 'move_extend', 'move_reduce', 'move_split',
            'move_merge', 'move_update', 'move_transfer', 'move_cross_over','move_identity']))
        for name in self.move_names.values():
            self.param_dict[name] = []
        self.track_samples = list()
        self.param_samples = defaultdict(list)
        self._log_pdf = list()
        self._acceptance_rate = list()
        # other control paramters
        self.do_record_samples = False
        self._tick = 0
        self.min_track_length = 0
        self.start_time = start_time

    def init_move_hist(self, move_histogram):
        '''This should be called just before starting to record samples'''
        self.last_hist = move_histogram

    def init_sample_stats(self, tick_time, sample_count):
        '''initial values to calc sample rate'''
        self.sample_rates = list()
        self.sample_counts = list()
        self._last_tick_time = tick_time
        self._last_sample_count = sample_count

    def set_min_track_length(self, min_track_length):
        '''for pair recording '''
        assert(isinstance(min_track_length, int))
        self.min_track_length = min_track_length

    def set_record_samples(self, true_or_false):
        self.do_record_samples = true_or_false

    def _record_parameters(self, state):
        ''' saving parameters '''
        params = json.loads(state.current_model_parameters.to_json())
        for k, v in params.items():
            if k not in self.param_dict:
                self.param_dict[k] = []
            self.param_dict[k].append(v)

    def _record_move_hist(self, state):
        '''saving the histogram. self.last_hist must be initialised'''
        move_hist = array(state.move_histogram)
        diff_hist, self.last_hist = move_hist - self.last_hist, move_hist
        for i in range(len(self.move_names)):
            self.param_dict[self.move_names[i]].append(diff_hist[i])

    def _record_sample(self, state):
        '''record partition sample'''
        indexed_partition = IndexedPartition()
        indexed_partition.from_partition(state.current_partition)
        self.track_samples.append([
            [track['observations'], track['time_span']] for track in indexed_partition.tracks()
            ])

    def _record_pairs(self, state):
        '''record averaged and non-averaged track pairs'''
        numtr = 0
        trlen = 0.0
        index_form = partition_to_index_form(state.current_partition)
        for track, idx_track in izip(state.current_partition.tracks, index_form):
            numtr += 1
            trlen += float(track.duration)
            vels_data = extract_kalman_states(track, state.current_model_parameters)
            vels_data = (vels_data, [0, 0, 0, 1], 1.0, [0, 0, 0, 1])

            if len(track.observations) < self.min_track_length:
                continue
            update_idx_pair_dicts(
                self.avg_pair_dict, self.pair_dict, track.observations, idx_track[0], vels_data, self._tick)
        if numtr > 0:
            self.param_dict['mean_track_length'].append(float(trlen)/float(numtr))
        else:
            self.param_dict['mean_track_length'].append(0.0)
        self.param_dict['track_count'].append(float(numtr))

    def _record_sample_stats(self, state):
        '''sample count and sample rate'''
        tick_time = time.time() - self.start_time
        sample_count = state.sample_count - self._last_sample_count
        sample_rate = float(sample_count)/float(tick_time - self._last_tick_time)
        self.sample_rates.append(sample_rate)
        self.sample_counts.append(sample_count)
        self._last_tick_time = tick_time
        self._last_sample_count = state.sample_count

    def record_state(self, state):
        '''records the sampling state'''
        if state.sample_count == 0:
            return
        self._tick += 1
        self._record_sample_stats(state)
        self._record_parameters(state)
        self._record_move_hist(state)
        self._log_pdf.append(state.current_log_pdf)
        self._acceptance_rate.append(state.recent_acceptance_rate)
        if self.do_record_samples:
            self._record_sample(state)
        else:
            self._record_pairs(state)

    def sample_stats(self):
        return self.sample_rates, self.sample_counts
    def parameter_items(self):
        return self.param_dict.items()
    def track_samples(self):
        return self.track_samples
    def log_pdf(self):
        return self._log_pdf
    def acceptance_rate(self):
        return self._acceptance_rate
    def pair_dict_items(self):
        return self.pair_dict.items()
    def avg_pair_dict_items(self):
        return self.avg_pair_dict.items()

def get_tracking_dir(input_fname):
    base_dir = os.path.dirname(input_fname)
    if os.path.basename(base_dir) == "input":
        base_dir = os.path.normpath(os.path.join(base_dir, ".."))
    directory = os.path.join(base_dir, "tracking")
    if not os.path.isdir(directory) and not os.path.exists(directory):
        os.mkdir(directory)
    if not os.path.isdir(directory):
        directory = base_dir
    return directory

def get_output_filename(input_fname, directory, infix = 'tracked', insert_time = True):
    '''determine the output file name
    FIXME: function returns different path depending on "insert_time"
    '''
    if not os.path.isdir(directory):
        raise Exception('Not a directory: "%s"' % directory)
    output = os.path.basename(input_fname).split('.')
    ext = output[-1]
    if not infix is None and len(infix):
        if infix in output:
            output = output[0:output.index(infix)]
            output.append(ext)
        if output[-2] == 'input':
            output[-2] = infix
        else:
            output.insert(-1, infix)
    else:
        if output[-2] == 'input':
            del(output[-2])
    output = '.'.join(output)
    if insert_time:
        (fname, ext) = os.path.splitext(os.path.join(directory, output))
        return fname + "." + datetime.utcnow().strftime('%y%b%d.%H%M') + ext
    else:
        return output


