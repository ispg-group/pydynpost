#!/usr/bin/env python
import os
import importlib
import numpy as np
import src.globalattr as glAttr
import src.griddata as grddt 
#import src.propagate as propagate
import src.traj as traj
import src.leaf as leaf

class Branch(glAttr.globalClass):

    def __init__(self, glbl, geom, time = None):

        super().__init__(glbl)
        self.geom = geom
        if hasattr(self, 'nrRNGs'):
            self.leaves = leaf.getSimpleIterator(glbl, geom) 
        #elif self.dynMethod == 'aims':
        #    self.nTraj, self.trajs = traj.getSimpleIterator(
        #        glbl, geom, rng=None, time=time
        #    ) 
        self.leaves = []

    def getMetric(self, metric):

        if hasattr(self, 'nrRNGs'): 
            return 0.0

        metricModule = 'src.' + metric
        metricModule = importlib.import_module(metricModule)
        getMetric = getattr(metricModule, 'get' + metric[0].upper() + metric[1:])
        metric = getMetric(self, str(self.geom))
        return metric

    def populateData(self):
        if self.dynMethod == 'aims':
            tmpTrajs = []
            for _traj in self.trajs:
                _traj.getStateID()
                _traj.getPositions()
                _traj.getMomenta()
                _traj.getWidths()
                _traj.getPhase()
                if self.model == 'zero':
                    _traj.pruneData()
                _traj.getAmplitude()
                error = _traj.getAmplitudeDot()
                if error != False:
                    _traj.calcAmplitudeDot()
                tmpTrajs.append(_traj)

            self.trajs = tmpTrajs
            
        else: 
            self.getPositions()
            self.getMomenta()
            self.getPhase()

    def saveNewAmps(self):
        for traj_ in self.trajs: 
            traj_.saveAmp()

def getSimpleIterator(glbl):

    def branchSimpleIterator():

        for geom in range(1, glbl.sampleSize + 1):
            if geom in glbl.dupList:
                continue

            currBranch = Branch(glbl, geom)
            yield currBranch
            
    simpleIterator = branchSimpleIterator()
         
    return simpleIterator 

def getPairIterator(glbl):
    
    def branchPairIterator(time = None):
        for geom1 in range(1, glbl.sampleSize + 1):
            if geom1 in glbl.dupList:
                continue

            for geom2 in range(geom1, glbl.sampleSize + 1):
                currBranch1 = Branch(glbl, geom1, time=time)
                currBranch2 = Branch(glbl, geom2, time=time)
                yield currBranch1, currBranch2
            
    pairIterator = branchPairIterator
         
    return pairIterator 

def fillGridData(branch):
    outPGrid = []
    outPot = []
    outCGrid = {}
    outCoup = {}

    for i in range(1,branch.nrStates):
        for j in range(i+1, branch.nrStates+1):
            outCGrid[f"c{i}_{j}"] =  []
            outCoup[f"c{i}_{j}"] =  []

    if hasattr(branch, 'nrRNGs'): 
        for leaf in branch.leaves:
            gridData = leaf.fillGridData([],[],[], [])   
            outPGrid.append(branch.filterDuplicates[gridData[0],Grid])
            outPot.append(branch.filterDuplicates[gridData[1],Pot])
            outCGrid.append(branch.filterDuplicates[gridData[2],Grid])
            outCoup.append(branch.filterDuplicates[gridData[3],Coup])

    else:
        gridData = getGridData(branch)
        outPGrid.extend(gridData[0])  
        outPot.extend(gridData[1])  
        for i in range(1,branch.nrStates):
            for j in range(i+1, branch.nrStates+1):
                outCGrid[f"c{i}_{j}"].extend(gridData[2][f"c{i}_{j}"])  
                outCoup[f"c{i}_{j}"].extend(gridData[3][f"c{i}_{j}"])  

    return (outPGrid, outPot, outCGrid, outCoup)

def getGridData(branch):
    return grddt.getGridData(branch.glbl, str(branch.geom))
