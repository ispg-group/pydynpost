#!/usr/bin/env python
import os
import importlib
import numpy as np
import src.globalattr as glAttr
import src.griddata as grddt 
import src.propagate as propagate
import src.leaf as leaf

class Branch(glAttr.globalClass):

    def __init__(self, glbl, geom):

        super().__init__(glbl)
        self.geom = geom
        if hasattr(self, 'nrRNGs'):
            self.leaves = leaf.getSimpleIterator(glbl, geom) 
        else:
            self.trajs = traj.getSimpleIterator(glbl, geom) 
        self.leaves = []

    def getMetric(self, metric):

        if hasattr(self, 'nrRNGs'): 
            return 0.0

        metric = 'src.' + metric
        metricModule = importlib.import_module(metric)
        getMetric = getattr(metricModule, 'get' + metric[0].upper() + metric[1:])
        metric = getMetric(self, str(self.geom))
        return metric

    def populateData(self):
        tmpTrajs = []
        for _traj in self.trajs:
            _traj.getID()
            _traj.getPositions()
            _traj.getMomenta()
            _traj.getPhase()
            error = _traj.getAmplitude()
            error = _traj.getDotAmplitude()
            tmpTrajs.append(_traj)
        
        self.trajs = tmpTrajs
        
        if error == True:
            self.postPropagate()

    def postPropagate(self):
        propagate(self) 
        self.saveNewAmps()

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
    
    def branchPairIterator():
        for geom1 in range(1, glbl.sampleSize + 1):
            if geom1 in glbl.dupList:
                continue

            for geom2 in range(geom1, glbl.sampleSize + 1):
                currBranch1 = Branch(glbl, geom1)
                currBranch2 = Branch(glbl, geom2)
                yield currBranch1, currBranch2
            
    pairIterator = branchPairIterator()
         
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
