# -*- coding: utf-8 -*-
"""
Ripley's K-function and L-function for observations

Created on Tue May  3 15:52:13 2016

@author: vcn81216
"""

import numpy as np
from scipy.spatial import distance

MAXSIZE = 6000

class SpaceFromTime(object):
    def __init__(self, pixel_per_frame):
        self.ppf = pixel_per_frame
        
    def __call__(self, frames):
        return self.ppf * frames
        
def ripley_k(observations, pixel_per_frame, d, area = None):
    '''calculates Ripley's k-function. Needs a (diffusion) rate in [pixel/frame] to convert time 
    into distance. If the area is *None* it will be calculated from the observations
    '''
    num_obs = len(observations)
    if num_obs > MAXSIZE:
        return dict()
    time_len = SpaceFromTime(pixel_per_frame)
    if isinstance(d, (int, float, long)):
        d = [d]
    if area is None:
        min_x, min_y, min_t = np.min(observations, axis = 0)
        max_x, max_y, max_t = np.max(observations, axis = 0)
        area = (max_x-min_x)*(max_y-min_y)*time_len(max_t-min_t)
    spaced_obs = [ [x, y, time_len(t)] for x, y, t in observations ]
    result = dict()
    for dist in d:
        cnt = 2 * sum(distance.pdist(spaced_obs)<dist)
        result[dist] = area /(num_obs*(num_obs-1)) * cnt
    return result
    
def ripley_l(observations, pixel_per_frame, d, area = None):
    '''Ripley's L-function'''
    return dict([[key, np.sqrt(val)/np.pi] for key, val in ripley_k(
        observations, pixel_per_frame, d, area).items()])
