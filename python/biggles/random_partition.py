# -*- coding: utf-8 -*-
"""
Create random partitions

Created on Thu May 19 11:47:48 2016

@author: vcn81216
"""

from _utilities import obs_dist
import numpy as np

SPEED_OF_LIGHT = 3.0

def valid_link(o1, o2):
    '''Is the link between two observations valid?
    The spatial distance must be <= SPEED_OF_LIGHT*(time difference)
    '''
    return obs_dist(o1, o2) <= SPEED_OF_LIGHT*(abs(o1[2]-o2[2])) 
    
def rev_partition(observations):
    '''create a random partition which is maximal using a greedy algorithm
    Start from the end
    '''
    obs_idx = range(len(observations))
    _, _, tmin = np.min(np.array(observations), axis = 0)
    _, _, tmax = np.max(np.array(observations), axis = 0)
    obs_dict = dict()
    tmin, tmax = int(tmin), int(tmax)
    for t in range(tmin, tmax +1):
        obs_dict[t] = []
    for idx in obs_idx:
        obs_dict[observations[idx][2]].append(idx)
    tracks = []
    for time in reversed(sorted(obs_dict.keys())):
        obs_at_t = list(obs_dict[time])
        np.random.shuffle(obs_at_t)
        for obs in obs_at_t:
            pre_track = [obs]
            done = False
            while not done:
                o = pre_track[-1]
                for dt in range(1, 5):
                    current_time = observations[o][2]
                    if current_time - dt < tmin:
                        break
                    cand = list(obs_dict[current_time - dt])
                    np.random.shuffle(cand)
                    for c in cand:
                        if valid_link(observations[o], observations[c]): # valid obs found
                            pre_track.append(c)
                            break
                    if o != pre_track[-1]: # there was an obs attached
                        break
                if o == pre_track[-1]: # no obs attached 
                    done = True
            if len(pre_track)>1: # i.e. track is valid
                pre_track.reverse()
                for o in pre_track:
                    obs_idx.remove(o)
                    obs_dict[observations[o][2]].remove(o)
                tracks.append({
                    'time_span' : [observations[pre_track[0]][2], observations[pre_track[-1]][2]+1],
                    'observations' : pre_track})
    return tracks, obs_idx
        
def max_partition(observations):
    '''create a random partition which is maximal using a greedy algorithm'''
    obs_idx = range(len(observations))
    _, _, tmin = np.min(np.array(observations), axis = 0)
    _, _, tmax = np.max(np.array(observations), axis = 0)
    obs_dict = dict()
    tmin, tmax = int(tmin), int(tmax)
    for t in range(tmin, tmax +1):
        obs_dict[t] = []
    for idx in obs_idx:
        obs_dict[observations[idx][2]].append(idx)
    tracks = []
    for time in sorted(obs_dict.keys()):
        obs_at_t = list(obs_dict[time])
        np.random.shuffle(obs_at_t)
        for obs in obs_at_t:
            pre_track = [obs]
            done = False
            while not done:
                o = pre_track[-1]
                for dt in range(1, 5):
                    current_time = observations[o][2]
                    if current_time + dt > tmax:
                        break
                    cand = list(obs_dict[current_time + dt])
                    np.random.shuffle(cand)
                    for c in cand:
                        if valid_link(observations[o], observations[c]): # valid obs found
                            pre_track.append(c)
                            break
                    if o != pre_track[-1]: # there was an obs attached
                        break
                if o == pre_track[-1]: # no obs attached 
                    done = True
            if len(pre_track)>1: # i.e. track is valid
                for o in pre_track:
                    obs_idx.remove(o)
                    obs_dict[observations[o][2]].remove(o)
                tracks.append({
                    'time_span' : [observations[pre_track[0]][2], observations[pre_track[-1]][2]+1],
                    'observations' : pre_track})
    return tracks, obs_idx
    
if __name__ == '__main__':
    import json
    import argparse
    parser = argparse.ArgumentParser(description=globals()['__doc__'])
    parser.add_argument('input', metavar='INPUT', type=str, 
        help='read partition from json file named INPUT')
    args = parser.parse_args()
    with open(args.input) as fh:
        data = json.load(fh)
    tracks, clutter = max_partition(data['observations'])
    total = len(data['observations'])
    oset = set(clutter)
    sub = len(clutter)
    
    for tr in tracks:
        oset |= set(tr['observations'])
        print("tr = %d" % len(tr['observations']))
    assert(len(oset) == total)
    print("tracks = {0}, clutter = {1}, total = {2}".format(
        len(tracks), len(clutter), total))
    