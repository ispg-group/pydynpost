#!/usr/bin/env python
import os
import importlib
import numpy as np
import globalattr as glAttr
import leaf

class Branch(glAttr.globalClass):

    def __init__(self, glbl, geom):

        super().__init__(glbl)
        self.geom = geom
        if hasattr(self, 'nrRNGs'):
            self.leaves = leaf.getSimpleIterator(glbl, geom) 
        self.glbl = glbl
        self.leaves = []

    def getMetric(self, metric):

        if hasattr(self, 'nrRNGs'): 
            return 0.0

        metricModule = importlib.import_module(metric)
        getMetric = getattr(metricModule, 'get' + metric[0].upper() + metric[1:])
        metric = getMetric(self.glbl, str(self.geom))
        return metric


def getSimpleIterator(glbl):

    def branchSimpleIterator():

        for geom in range(1, glbl.sampleSize + 1):
            if geom in glbl.dupList:
                continue

            currBranch = Branch(glbl, geom)
            yield currBranch
            
    simpleIterator = branchSimpleIterator()
         
    return simpleIterator 
