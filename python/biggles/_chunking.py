# -*- coding: utf-8 -*-
"""
Chunking utilities for Biggles

Created on Thu Mar 31 10:24:48 2016

@author: vcn81216
"""

import numpy as np

class GenericChunk(object):
    '''a chunk'''
    def __init__(self, loc, t, path, apron, size):
        self.loc = loc
        self.t = t
        self.path = path
        self._apron = apron
        self._size = size

    def __len__(self):
        return len(self.data)
    def __str__(self):
        return 'DataChunk [[{0:3d}, {1:3d}], [{2:3d}, {3:3d}], [{4:.0f}, {5:.0f}]], count = {6:6d}. path = "{7}"'.format(
            int(self.loc[0][0]), int(self.loc[0][1]), 
            int(self.loc[1][0]), int(self.loc[1][1]), 
            self.t[0], self.t[1],
            len(self.data), self.path)
        
    def extent(self, ax):
        return self.loc[ax][1]-self.loc[ax][0]
        
    def area(self):
        return self.extent(0) * self.extent(1)

    def set_apron(self, new_apron_size):
        assert(new_apron_size>=0)
        self._apron = new_apron_size
        
    def apron(self):
        return self._apron
        
    def set_size(self, new_size):
        assert(new_size>0)
        self._size = new_size
        
    def size(self):
        return self._size
        

class DataChunk(GenericChunk):
    '''a chunk of (clutter) observations'''
    def __init__(self, data, loc, t, path, apron = 0, size = set("over write me")):
        super(DataChunk, self).__init__(loc, t, path, apron, size)
        self.data = data
        
    def cropped(self, loc):
        '''return spatially cropped data'''
        rmin = np.array([loc[0][0], loc[1][0], self.t[0]])
        rmax = np.array([loc[0][1], loc[1][1], self.t[1]])
        return [ obs for obs in self.data if (rmin <= obs).all() and (obs < rmax).all() ]
        
    def get_chunks(self, chunklist, ax = 0):
        ''' Splits *chunk* into small chunk of max size *self.size()* and appends them to chunklist.
        This will cut the chunk in half as long as it contains more observations than *self.size()*.
        The cut alternately goes along the x-axis and y-axis
        '''
        if len(self) <= self.size() or self.extent(ax) <= 2.*self._apron:
            chunklist.append(self)
            #print(self)
            return
        locmed = (self.loc[ax][1]+self.loc[ax][0])/2.0
        loclo = list(self.loc)
        loclo[ax] = [self.loc[ax][0], locmed + self._apron/2.0]
        lochi = list(self.loc)
        lochi[ax] = [locmed - self._apron/2.0, self.loc[ax][1]]
        t = self.t
        DataChunk(self.cropped(loclo), loclo, t, 
            self.path + "lb"[ax], self.apron(), self.size()).get_chunks(chunklist, (ax+1)%2)
        DataChunk(self.cropped(lochi), lochi, t, 
            self.path + "rt"[ax], self.apron(), self.size()).get_chunks(chunklist, (ax+1)%2)
        
    @staticmethod
    def from_partition(partition):
        '''creates a Data chunk from a partition'''
        data = partition['clutter']
        for track in partition['tracks']:
            data.extend(track['observations'])
        xmin, ymin, tmin = np.floor(np.min(data, axis=0)).astype(float)
        xmax, ymax, tmax = np.ceil(np.max(data, axis=0)).astype(float) + 1
        return DataChunk(data, [[xmin, xmax], [ymin, ymax]], [tmin, tmax], "")
           
class IndexedDataChunk(GenericChunk):
    '''a chunk of an index partition'''
    def __init__(self, observations, data, loc, t, path, apron = 0, size = set("over write me")):
        super(IndexedDataChunk, self).__init__(loc, t, path, apron, size)
        self.obs_list = observations # the observations [[x, y, t], ...]
        self.data = data # observation indices
        
    def _ge(self, index, obs):
        '''is the observation given by *index* greater or equal the numpy array *obs*?'''
        return (obs <= self.obs_list[index]).all()
    def _lt(self, index, obs):
        '''is the observation given by *index* less than the numpy array *obs*?'''
        return (obs > self.obs_list[index]).all()
    def cropped(self, loc):
        '''return spatially cropped data'''
        rmin = np.array([loc[0][0], loc[1][0], self.t[0]])
        rmax = np.array([loc[0][1], loc[1][1], self.t[1]])
        return [ idx for idx in self.data if self._ge(idx, rmin) and self._lt(idx, rmax) ]
        
    def get_chunks(self, chunklist, ax = 0):
        ''' Splits *chunk* into small chunk of max size *size* and appends them to chunklist.
            This will cut the chunk in half as long as it contains more observations than *size*.
            The cut alternately goes along the x-axis and y-axis
        '''
        if len(self) <= self.size() or self.extent(ax) <= 2.*self._apron:
            chunklist.append(self)
            #print(self)
            return
        locmed = (self.loc[ax][1]+self.loc[ax][0])/2.0
        loclo = list(self.loc)
        loclo[ax] = [self.loc[ax][0], locmed + self._apron/2.0]
        lochi = list(self.loc)
        lochi[ax] = [locmed - self._apron/2.0, self.loc[ax][1]]
        t = self.t
        IndexedDataChunk(self.obs_list, self.cropped(loclo), loclo, t, 
            self.path + "lb"[ax], self.apron(), self.size()).get_chunks(chunklist, (ax+1)%2)
        IndexedDataChunk(self.obs_list, self.cropped(lochi), lochi, t, 
            self.path + "rt"[ax], self.apron(), self.size()).get_chunks(chunklist, (ax+1)%2)
        
    @staticmethod
    def from_partition(partition):
        '''creates a Data chunk from a indexed partition'''
        obs_list = np.array(partition['observations'])
        data = partition['clutter']
        for track in partition['tracks']:
            data.extend(track['observations'])
        xmin, ymin, tmin = np.floor(np.min(obs_list, axis=0)).astype(float)
        xmax, ymax, tmax = np.ceil(np.max(obs_list, axis=0)).astype(float) + 1
        return IndexedDataChunk(obs_list, data, [[xmin, xmax], [ymin, ymax]], [tmin, tmax], "")
