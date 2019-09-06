# -*- coding: utf-8 -*-
"""
hdf5 interface for biggles

Created on Tue May 10 15:04:38 2016

@author: vcn81216
"""

import json
from numpy import array
from biggles import IndexedPartition

def hdf5tracks3(h5tracks):
    '''reverses  track3hdf5'''
    assert('enc' in h5tracks.attrs)
    assert(h5tracks.attrs['enc'] == u"ttlo..")
    tracks = h5tracks[:].tolist()
    result = []
    while len(tracks):
        nobs = tracks[2]
        track = {'time_span' : tracks[0:2], 'observations' : tracks[3:3+nobs] }
        tracks = tracks[3+nobs:]
        result.append(track)
    return result
            
def tracks2hdf5(tracks, root):
    '''write tracks to hdf5'''
    for i, track in enumerate(tracks):
        h5tr = root.create_dataset('{0:d}'.format(i), data = array(
            track['observations'], dtype = float))
        h5tr.attrs['dur'] = array(track['time_span'], dtype = int)
        
def tracks3hdf5(tracks, root, name = 'tracks'):
    '''Transforms indexed tracks and writes them to hdf5.
    The tracks will be transformed into a single data set. 1) 2 ints for the time span
    2) 1 int, N,  for the number of observations 3) N ints for the observations 
    '''
    trans = list()
    for track in tracks:
        trans.extend(track['time_span'])
        trans.append(len(track['observations']))
        trans.extend(track['observations'])
    dset = root.create_dataset(name, data = array(trans, dtype = int))
    dset.attrs['enc'] = u"ttlo.."
        
def truth2hdf5(truth, root):
    '''writes ground truth to hdf5'''
    for i, gt in enumerate(truth):
        h5gt = root.create_dataset('%d' % i, data = array(gt[1], dtype = float) )
        h5gt.attrs['t0'] = gt[0]

def write_hdf5_partition(partition, root):
    '''writes a partition to a hdf5 file''' 
    part = json.loads(partition.to_json())
    root.create_dataset('clutter', data = array(part['clutter'], dtype = float))
    tracks2hdf5(part['tracks'], root.create_group("tracks"))
    
def write_hdf5_observations(observations, root):
    '''writes observations to a hdf5 file'''
    root.create_dataset('observations', data = array(observations, dtype = float))
    
def write_hdf5_indexed_partition(part, root):
    '''writes a partition to a hdf5 file'''
    assert(isinstance(part, IndexedPartition)) # just so that you know
    root.create_dataset('clutter', data = array(part.clutter(), dtype = int))
    tracks3hdf5(part.tracks(), root.create_group("tracks"))
    
def write_hdf5_parameters(para, root):
    '''writes the biggles model parameters'''
    root.attrs['births_per_frame'] = para.mean_new_tracks_per_frame
    root.attrs['clutter_per_frame'] = para.mean_false_observations_per_frame
    root.attrs['observation_probability'] = para.generate_observation_probability
    root.attrs['survival_probability'] = para.frame_to_frame_survival_probability
    oec = array(para.observation_error_covariance, dtype = float)
    root.attrs['observation_covariance'] = oec

def write_hdf5_track_samples(track_samples, root):
    '''writes track samples to hdf5
    the input data consists of a list of samples, each sample is a partition: a list of tracks
    each track is a list with tow elements: a list of observation and a two-element-list 
    containing the time span. This is somewhat a doubling of "tracks3hdf5"
    '''
    for i, part in enumerate(track_samples):
        trans = list()
        for track in part:
            trans.extend(track[1]) # 2 ints - time span
            trans.append(len(track[0])) # 1 int - num obs
            trans.extend(track[0]) # some ints - obs
        dset = root.create_dataset(str(i), data = array(trans, dtype = int))
        dset.attrs['enc'] = u"ttlo.." # encoding descriptor        
    
def write_hdf5_samples(samples, root):
    '''writes biggles track samples to hdf5'''
    root.create_dataset('log_pdf', data = array(samples['log_pdf'], dtype = float))
    root.attrs['sample_count'] = samples['sample_count']
    root.attrs['burnin_length'] = samples['burnin_length']
    root.attrs['chunk_id'] = samples['chunk_id']
    root.attrs['grid_row'] = samples['grid_row']
    root.attrs['record_ticks'] = samples['record_ticks']
    root.attrs['x_extent'] = array(samples['x_extent'], dtype = float)
    root.attrs['y_extent'] = array(samples['y_extent'], dtype = float)
    root.create_dataset('acceptance_rate', data = array(samples['acceptance_rate'], dtype = float))
    h5para = root.create_group('param_samples')
    for key, value in samples['param_samples'].items():
        h5para.create_dataset(key, data = array(value, dtype = float))
    write_hdf5_track_samples(samples['track_samples'], root.create_group('track_samples'))

def write_hdf5_sim(data, fh):
    '''writes data to a hdf5 file'''
    fh.attrs['simulation'] = json.dumps(data['simulation'])
    fh.attrs['fileheader'] = json.dumps(data['fileheader'])
    fh.attrs['metadata'] = json.dumps(data['metadata'])
    fh.attrs['flags'] = json.dumps(data['flags'])
    fh.create_dataset('observations', data = array(data['observations'], dtype = float))
    fh.create_dataset('clutter', data = array(data['clutter'], dtype = int))
    tracks3hdf5(data['tracks'], fh)
    truth2hdf5(data['ground_truth'], fh.create_group('ground_truth'))
    truth2hdf5(data['unobserved'], fh.create_group('unobserved'))
