# -*- coding: utf-8 -*-
"""
general utilities

Created on Fri Apr  1 08:56:28 2016

@author: vcn81216
"""

import numpy as np
from collections import deque
import os

def obs_dist(o1, o2):
    '''distance between observations at the same time point'''
    return np.sqrt((o1[0] - o2[0])**2 + (o1[1] - o2[1])**2)

def sec2str(seconds):
    '''converts seconds to "h:mm:ss"'''
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    return "%d:%02d:%02d" % (h, m, s)

def count2str(count):
    '''converts a count value into a string'''
    if count < 1000:
        return "{0:6d}".format(count)
    index = 0
    num = float(count)
    while num >= 999.95:
        num /= 1000.0
        index += 1
    return "{0:5.1f}{1}".format(num, ' kMGTP'[index])

def four2str(count):
    '''converts a count value into a 4 character string'''
    if count < 9.5:
        return "{0:4.2f}".format(count)
    if count < 99.5:
        return "{0:4.1f}".format(count)
    if count < 10000:
        return "{0:4.0f}".format(count)
    index = 0
    num = float(count)
    while num >= 999.5:
        num /= 1000.0
        index += 1
    return "{0:3.0f}{1}".format(num, ' kMGTP'[index])

def float2str(value):
    '''converst a float value into a string'''
    if abs(value) < 1.0:
        index = 0
        while abs(value) < 1.0:
            value *= 1000.0
            index += 1
        formstr = u"{0:6.%df}{1}" % (4-int(np.floor(np.log10(abs(value))+1)))
        return formstr.format(value, u" m\u03bcnpfa"[index])
    if abs(value) < 1000:
        return "{0:7.1f}".format(value)
    index = 0
    while abs(value) >= 1000.0:
        value /= 1000.0
        index += 1
    formstr = "{0:6.%df}{1}" % (4-int(np.floor(np.log10(abs(value))+1)))
    return formstr.format(value, ' kMGTPE'[index])

def percentstr(ratio):
    '''takes a ratio returns a 3 char string repr the percentage'''
    assert(0.0 <= ratio < 10)
    tmp = 100.0*ratio
    if tmp >= 9.5:
        return "{0:3.0f}".format(tmp)
    else:
        return "{0:3.1f}".format(tmp)

def autocorr(data, i):
    ''' auto correlation of *data* with lag *i*
    '''
    return np.corrcoef(data[i:], data[:-i])[1,0]

def observation_count(partition):
    return len(partition.clutter) + sum([len(t.observations) for t in partition.tracks])

def obscmp(o1, o2):
    '''compares two observations given as list of 3'''
    if o1[2] < o2[2]: return -1
    if o1[2] > o2[2]: return 1
    if o1[0] < o2[0]: return -1
    if o1[0] > o2[0]: return 1
    if o1[1] < o2[1]: return -1
    if o1[1] > o2[1]: return 1
    return 0

def get_extent(partition):
    obs = list(partition['clutter'])
    for t in partition['tracks']:
        obs.extend(t.observations)
    obs = np.array(obs, dtype = float)
    xmin, ymin, tmin = np.floor(np.min(obs, axis=0))
    xmax, ymax, tmax = np.ceil(np.max(obs, axis=0))
    return [xmin, xmax], [ymin, ymax], [tmin, tmax], len(obs)

def fifo(length = 2000):
    return deque(maxlen = length)

def number_parser(num_str, num_type=float):
    '''converts number strings using metric prefixes into numbers'''
    postfixes = " kMGTPE"
    power = 0
    if num_str[-1] in postfixes:
        power = postfixes.index(num_str[-1])
        num_str = num_str[:-1]
    return num_type(float(num_str.replace(',','')) * 1000**power)

def auto_corr_fun(y):
    if len(y) < 2:
        return list()
    yy = np.array(y)
    yyunbiased = yy - np.mean(yy)
    yynorm = np.sum(yyunbiased**2)
    if yynorm == 0.0:
        return list()
    acor = np.correlate(yyunbiased, yyunbiased, 'same')/yynorm
    return acor[len(acor)/2:]

def sample_pairs(size_data, num_samples):
    '''uniformly samples "num_samples" unique pairs (x, y) so that 0 <= x < y < "size_data"'''
    assert(size_data > 0)
    assert(0 <= num_samples < (size_data -1)*size_data/2)
    result = set()
    pool = range(size_data)
    while len(result) < num_samples:
        cand = np.random.choice(pool, 2, replace = False)
        result.add((min(cand), max(cand)))
    return result

def auto_corr_time(y):
    acf = auto_corr_fun(y)
    tot = 0
    for rho in acf[1:]:
        if rho < 0.1:
            break
        tot += rho
    return 1.0 + 2.0*tot

def get_random_string(num):
    return ''.join(np.random.choice(list("abcdefghijklmnopqrstuvwxyz0123456789"), num).tolist())

class UnsafeOpenError(Exception):
    pass

def safe_open(fname, mode='r', unsafe=False, logger = lambda x : None):
    if not unsafe:
        if os.path.exists(fname):
            logger('%s: will not overwrite existing file unless --unsafe option is specified.' % (fname,))
            raise UnsafeOpenError()
    return open(fname, mode)

