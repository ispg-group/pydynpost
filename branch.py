#!/usr/bin/env python
import os
import importlib
import numpy as np
import globalattr as glAttr

class Branch(glAttr.globalClass):
    def __init__(self, glbl, geom, rng=None):
        super().__init__(glbl)
        self.geom = geom
        self.rng = rng

    def getMetric(self, metric):
        metricModule = importlib.import_module(metric)
        getMetric = getattr(metricModule, 'get' + metric)
        if self.rng != None:
            metric = getMetric(self.dynMethod, self.pckg, self.CWD,
                               self.geom, self.geomDir, self.rng,
                               self.rngDir)
        else:
            metric = getMetric(self.dynMethod, self.pckg, self.CWD,
                               self.geom, self.geomDir)
        return metric
        

def getSimpleIterator(glbl):
    rngLoop = False
    if glbl.dynMethod == 'tsh': 
        if hasattr(glbl, 'nrRNGs'):
            rngLoop = True
    elif glbl.dynMethod == 'aims':
        if glbl.AIMStype != 'AIMS':
            assert hasattr(glbl, 'nrRNGs') 
            rngLoop = True

    def branchSimpleIterator():
        for geom in range(1, glbl.sampleSize + 1):
            if geom in glbl.dupList:
                continue

            if rngLoop:
                for rng in range(1, glbl.nrRNGs + 1):
                    print(geom, rng)
                    currBranch = Branch(glbl, geom, rng=rng)
                    yield currBranch 
            else:
                currBranch = Branch(glbl, geom)
                yield currBranch 
            
    bIterator = branchSimpleIterator()
         
    return bIterator 


