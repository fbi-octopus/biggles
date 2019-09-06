# -*- coding: utf-8 -*-
"""
Partition vector for enabling partiton arithmetic

Created on Thu Sep 22 15:12:58 2016 @author: vcn81216
"""

from collections import defaultdict
from scipy.sparse import dok_matrix, isspmatrix

class PV2(object):
    def __init__(self, init):
        if isinstance(init, (int, long)) and init>0:
            self._data = dok_matrix((init, init))
        elif isinstance(init, PV2):
            self._data = init._data
        elif isspmatrix(init):
            self._data = dok_matrix(init)
        else:
            raise TypeError("unsupported operand")
    @property
    def data(self):
        return self._data
    def from_indexed_partition(self, indexed_partition):
        for track in indexed_partition.tracks():
            for o1, o2 in zip(track['observations'][0:-1], track['observations'][1:]):
                self._data[o1, o2] += 1
    def __add__(self, other):
        return PV2(self.data + other.data)
    def __sub__(self, other):
        return PV2(self.data - other.data)
    def _dot_product(self, other):
        return self.data.multiply(other.data).sum()
        
    def _scale(self, other):
        return PV2(self.data*other)

    def __mul__(self, other):
        if isinstance(other, PV2):
            return self._dot_product(other)
        if isinstance(other, (float, int, long)):
            return self._scale(other)
        raise TypeError("unsupported operand")
        
    def __div__(self, other):
        return self.__truediv__(other)
        
    def __truediv__(self, other):
        return PV2(self.data/other)
        
    def __pow__(self, other, modulo = None):
        return PV2(self.data**other)

class PartitionVector(object):
    def __init__(self, init = None, type = float):
        assert(type in [float, int, long])
        self._type = type
        self.clear()
        if isinstance(init, PartitionVector):
            self._type = init._type
            self._data = dict(init._data)
        elif not init is None:
            raise TypeError("unsupported operand")

    def clear(self):
        self._data = defaultdict(self._type)
        
    def keys(self):
        return self._data.keys()
    
    def from_indexed_partition(self, indexed_partition):
        self.clear()
        for track in indexed_partition.tracks():
            for o1, o2 in zip(track['observations'][0:-1], track['observations'][1:]):
                self[(o1, o2)] += 1
                
    def __getitem__(self, key):
        return self._data[key]
    def __setitem__(self, key, value):
        self._data[key] = value
    
    def __add__(self, other):
        result = PartitionVector()
        keys = set(self.keys()) | set(other.keys())
        for key in keys:
            result[key] = self[key] + other[key]
        return result

    def __sub__(self, other):
        result = PartitionVector()
        keys = set(self.keys()) | set(other.keys())
        for key in keys:
            result[key] = self[key] - other[key]
        return result
        
    def _dot_product(self, other):
        return sum([
            self[key]*other[key] for key in set(self.keys()) | set(other.keys())
        ])
        
    def _scale(self, other):
        result = PartitionVector(self)
        for key in self.keys():
            result[key] *= other
        return result

    def __mul__(self, other):
        if isinstance(other, PartitionVector):
            return self._dot_product(other)
        if isinstance(other, (float, int, long)):
            return self._scale(other)
        raise TypeError("unsupported operand")
        
    def __div__(self, other):
        return self.__truediv__(other)
        
    def __truediv__(self, other):
        if isinstance(other, (float, int, long)):
            result = PartitionVector(self)
            for key in self.keys():
                result[key] /= other
            return result
        raise TypeError("unsupported operand")
        
    def __pow__(self, other, modulo = None):
        if other == 2 :
            return self._dot_product(self)
        raise TypeError("unsupported operand")

