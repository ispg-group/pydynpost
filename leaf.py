#!/usr/bin/env python
import os
import importlib
import numpy as np
import globalattr as glAttr

class Leaf(glAttr.globalClass):
    def __init__(self, glbl, geom, rng):
        super().__init__(glbl)
        self.geom = geom
        self.geom = rng

    def getMetric(self, metric):
        metricModule = importlib.import_module(metric)
        getMetric = getattr(metricModule, 'get' + metric)
        metric = getMetric(self.dynMethod, self.pckg, self.CWD,
                           self.geom, self.geomDir, self.rng,
                           self.rngDir)
        return metric


def getSimpleIterator(glbl, geom):

    def leafSimpleIterator():
        for rng  in range(1, glbl.nrRNGs + 1):
            currLeaf = Leaf(glbl, geom, rng)
            print(rng, geom)
            yield currLeaf
    
    return leafSimpleIterator() 
