# -*- coding: utf-8 -*-
"""
interface to convert partitions to indexed partitions

Created on Tue Apr 19 14:58:05 2016

@author: vcn81216
"""

class Flag(object):
    INDEXED, OBSERVATIONS = range(2)

def indexed_data(tracks, clutter):
    '''converts clutter and tracks given as lists of observations 
    into list of indices 
    return the list observations and indexed tracks and clutter
    '''
    observations = list(clutter)
    indexed_clutter = range(len(clutter))
    indexed_tracks = list()
    for t in tracks :
        first = len(observations)
        observations.extend(t['observations'])
        indexed_tracks.append({
            'time_span' : t['time_span'],
            'observations': range(first, len(observations))
        })
    return observations, indexed_tracks, indexed_clutter
        
def observation_data(observations, indexed_tracks, index_clutter):
    '''converts indexed clutter and indexed tracks into lists obs of observations
    Returns (*tracks*, *clutter*)
    '''
    return [ { 'time_span' : t['time_span'], 
               'observations' : [   observations[i] for i in t['observations'] ]
           } for t in indexed_tracks ], [ observations[i] for i in index_clutter ]

def indexed_partition(partition):
    '''converts a observation partition into an index partition'''
    new_part = dict(partition)
    new_part['observations'], new_part['tracks'], new_part['clutter'] = indexed_data(
        partition['tracks'], partition['clutter']
    )
    replace_flag(new_part, Flag.OBSERVATIONS, Flag.INDEXED)
    return new_part
    
def has_flag(partition, flag):
    '''tests if a data set has a flag'''
    return 'flags' in partition and flag in partition['flags']
    
def replace_flag(partition, old_flag, new_flag):
    '''modify partition flags'''
    if not 'flags' in partition:
        partition['flags'] = []
    if old_flag in partition['flags']:
        partition['flags'].remove(old_flag)
    partition['flags'].append(new_flag)

def observation_partition(partition):
    '''converts an indexed partition into a observation partition'''
    new_partition = dict(partition)
    del new_partition['observations']
    new_partition['tracks'], new_partition['clutter'] = observation_data(
        partition['observations'], partition['tracks'], partition['clutter'] 
    )
    replace_flag(new_partition, Flag.INDEXED, Flag.OBSERVATIONS)
    return new_partition