# -*- coding: utf-8 -*-
"""
Tracking utilities

Created on Fri Apr 15 10:11:35 2016

@author: vcn81216
"""

from __future__ import print_function
from numpy import array, append, hypot, mean
from itertools import izip, islice
from biggles import ModelParameters, IndexedPartition, observation_pool
import json
import time
from collections import defaultdict

class Tracker(object):
    def __init__(self, stepper):
        self.stepper = stepper
        self.sample_rate = 0

    def run(self, steps):
        t0 = time.time()
        c0 = self.stepper.tracking_state.sample_count
        self.stepper.step(steps)
        self.sample_rate = float(self.stepper.tracking_state.sample_count - c0)/float(time.time() - t0)

        return self.state()

    def sample_observation_error(self, lag):
        self.stepper.sample_observation_error(lag)
    def fix_observation_error(self, R):
        self.stepper.fix_observation_error(R*R, 0, R*R)
    def state(self):
        return self.stepper.tracking_state

class TickTimer(object):
    def __init__(self, lag):
        self._lag = lag
        self._last = 0

    def passed(self, time):
        if time < self._last + self._lag:
            return False
        self._last = time
        return True

def update_square_diffs(track, res):
    ''' adds the square displacements to the dictionary *res*.
        track is a list of observations
    '''
    num_obs = len(track)
    for i in range(num_obs-1):
        o1 = track[i]
        for j in range(i+1, num_obs):
            o2 = track[j]
            res[o2[2] - o1[2]].append((o2[0] - o1[0])**2 + (o2[1] - o1[1])**2)

def msd_partition(tracks, cutoff = 50):
    '''the mean square displacement function *tracks*'''
    square_disp = defaultdict(list)
    for track in tracks:
        update_square_diffs(track, square_disp)
    msd = dict()
    for key, value in square_disp.items():
        if key <= cutoff:
            msd[key] = mean(value)
    return msd.items()

def states_to_tracks(state_tracks):
    ''' Formally converts tracks of states into tracks of observations. Its sole purpose
        is to avoid doubling functionality for observation tracks and state tracks
    '''
    tracks = []
    for t, strack in state_tracks:
        otrack = {
            'observations' : [],
            'time_span' : [t, t+len(strack)]
        }
        for state in strack:
            otrack['observations'].append([state[0], state[2], t])
            t += 1
        tracks.append(otrack)
    return tracks

def extract_kalman_states(track, current_model_parameters):
    ''' | Return the the state estimates for the observations of a track
        | vels = [
        |   1.0,
        |   state: x velocity, state: y velocity,
        |   state variance: x velocity, state variance: y velocity,
        |   state: x postion, state: y postion,
        |   state variance: x postion, state variance: y postion,
        | ]
    '''
    vels = []
    for sc in track.make_kalman_filter(current_model_parameters).corrections:
        s = sc.state # qx, vx, qy, vy
        c = sc.covariance
        ''' length of each vector: 9 '''
        vels.append(array([1, s[1], s[3], c[1][1], c[3][3], s[0], s[2], c[0][0], c[2][2]]))
    return vels


def update_pair_dicts(avg_pair_dict, pair_dict, obs, vels_data, record_tick):
    ''' updates the pair dictionaries
    '''


    vels, pbins, tstat, _ =  vels_data
    ''' average data:  9 + 1 + 2 + 3
        9 = 1 + 8
            1 always "1"
            8 kalman filter states + variances
        1 t-statistics for velocity (disused)
        2 number of observations & number of time stamps between observations
        3 p-value bins for t-test
    '''
    avg_vels = append( mean(array(vels), axis=0), [
        tstat,
        len(obs), obs.last_time_stamp - obs.first_time_stamp,
        pbins[0], pbins[1], pbins[2]
        ])
    assert avg_vels[0] == 1.0
    assert(len(avg_vels) == 15)

    for o1, o2, v1, v2 in izip( obs, islice(obs, 1, None), vels, islice(vels, 1, None) ):
        key = ((o1.x, o1.y, o1.t), (o2.x, o2.y, o2.t))
        ''' appending to not averaged:
            - the distance to the next state,
            - the distance estimated state to observation
            length now: 9+2 = 11
        '''
        v = append(v1, [ hypot(v1[5]-v2[5], v1[6]-v2[6])/(o2.t-o1.t),
                           hypot(v1[5]-o1.x,v1[6]-o1.y) ])
        try:
            pair_dict[key] = pair_dict[key] + v
        except KeyError:
            pair_dict[key] = v

        assert pair_dict[key][0] <= record_tick
        try:
            avg_pair_dict[key] = avg_pair_dict[key] + avg_vels
        except KeyError:
            avg_pair_dict[key] = avg_vels

        assert avg_pair_dict[key][0] <= record_tick

def update_idx_pair_dicts(avg_pair_dict, pair_dict, obs, obs_idx, vels_data, record_tick):
    ''' updates the pair dictionaries
    '''
    vels, pbins, tstat, _ =  vels_data
    ''' average data:  9 + 1 + 2 + 3
        9 = 1 + 8
            1 always "1"
            8 kalman filter states + variances
        1 t-statistics for velocity (disused)
        2 number of observations & number of time stamps between observations
        3 p-value bins for t-test
    '''
    avg_vels = append( mean(array(vels), axis=0), [
        tstat,
        len(obs), obs.last_time_stamp - obs.first_time_stamp,
        pbins[0], pbins[1], pbins[2]
        ])
    assert avg_vels[0] == 1.0
    assert(len(avg_vels) == 15)

    assert(len(obs) == len(obs_idx))

    idx_dict = dict(zip(obs_idx, obs))

    for io1, io2, v1, v2 in izip( obs_idx, islice(obs_idx, 1, None), vels, islice(vels, 1, None) ):
        #key = ((o1.x, o1.y, o1.t), (o2.x, o2.y, o2.t))
        key = (io1, io2)
        ''' appending to not averaged:
            - the distance to the next state,
            - the distance estimated state to observation
            length now: 9+2 = 11
        '''
        o1, o2 = idx_dict[io1], idx_dict[io2]
        v = append(v1, [ hypot(v1[5]-v2[5], v1[6]-v2[6])/(o2.t-o1.t),
                           hypot(v1[5]-o1.x,v1[6]-o1.y) ])
        try:
            pair_dict[key] = pair_dict[key] + v
        except KeyError:
            pair_dict[key] = v

        assert pair_dict[key][0] <= record_tick
        try:
            avg_pair_dict[key] = avg_pair_dict[key] + avg_vels
        except KeyError:
            avg_pair_dict[key] = avg_vels

        assert avg_pair_dict[key][0] <= record_tick

def observation_velocities(obs):
    ''' returns the velocities along a track by calculating the displacement of
        consecutive observations divided by the time lag.
        returns array of pairs of the velocities
        [ (vx1, vy1), (vx2, vy2), ...]
    '''
    return array([
            [ (float(o2.x-o1.x)/float(o2.t-o1.t)), (float(o2.y-o1.y)/float(o2.t-o1.t)) ]
            for o1, o2 in izip( obs, islice(obs, 1, None) )
        ])

def init_parameters():
    '''initialise the model parameters'''
    return ModelParameters.from_json(json.dumps({
            'births_per_frame': 0.05,
            'clutter_per_frame': 0.05,
            'survival_probability': 0.95,
            'observation_probability': 0.7,
            'observation_covariance': [[0.1, 0], [0, 0.1]]
    }))

def est_max_ged(observations):
    ''' Estimation of the maximum graph edit distance (1st approx)
    for each pair of consecutive frames add thd smaller number of obs to the estimate
    '''
    hist = defaultdict(int)
    for x, y, t in observations:
        hist[t] += 1
    times = hist.keys()
    times.sort()
    result = 0
    for fst, snd in zip(times[0:-1], times[1:]):
        result += min(hist[fst], hist[snd])
    return result

def get_max_ged(observations):
    ''' Maximum graph edit distance. No maximum spatial distance is assumed.
    A link is considered valid if the observations have different time stamps
    '''
    hist = defaultdict(int)
    for _, _, t in observations:
        hist[t] += 1
    times = hist.keys()
    times.sort()
    result = 0
    carried = 0
    for fst, snd in zip(times[0:-1], times[1:]):
        c1 = hist[fst] + carried
        c2 = hist[snd]
        result += min(c1, c2)
        carried = max(0, c1-c2)
    return result

def convert_initial_partition(initial_partition, pool = None):
    '''Converts a indexed partition in dict format into a Biggles partition.
    If the inital partition also contains the realisation of a simulated data
    it also returns a biggles partition of the realisation.
    Returns (bool "have simulation", Partition "initial partition", Partition "simulation",
             ObservationPool "observations")
    '''

    if pool is None:
        pool = observation_pool()
        for x, y, t in initial_partition['observations']:
            pool.push(float(x), float(y), int(t))
    else:
        assert(isinstance(pool, type(observation_pool())))
    partition_adapter = IndexedPartition()
    partition_adapter.set_clutter(initial_partition['clutter'])
    partition_adapter.set_tracks(initial_partition['tracks'])
    have_sims = False
    sims_partition = None
    if "initial_partition" in initial_partition:
        have_sims = True
        sims_adapter = IndexedPartition()
        sims_adapter.set_clutter(initial_partition['initial_partition']['clutter'])
        sims_adapter.set_tracks(initial_partition['initial_partition']['tracks'])
        sims_partition = sims_adapter.to_partition(pool)
    initial_partition = partition_adapter.to_partition(pool)
    return (have_sims, initial_partition, sims_partition, pool)

