#!/usr/bin/env python
import os
import importlib
import numpy as np
import globalattr as glAttr

class branch(glAttr.globalClass):
    def __init__(self, glbl, geom, rng=None):
        super().__init__(glbl)
        self.geom = geom
        if rng != None:
            self.rng = rng

    def getMetric(self, metric):
        print(metric, self.geom, self.rng)
        

def getSimpleIterator(glbl):
    rngLoop = False
    if glbl.dynMethod == 'tsh': 
        assert hasattr(glbl, 'nrRNGs') 
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
                    currBranch = branch(glbl, geom, rng=rng)
                    yield currBranch 
            else:
                currBranch = branch(glbl, geom)
                yield currBranch 
            
    bIterator = branchSimpleIterator()
         
    return bIterator 


