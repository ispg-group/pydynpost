#!/usr/bin/env python
import os
import importlib
import itertools
import numpy as np
import globalattr as glAttr
import branch

class Tree(glAttr.globalClass):

    def __init__(self, glbl):

        super().__init__(glbl)
        self.branches     = branch.getSimpleIterator(glbl)
        self.branchPairs  = branch.getPairIterator(glbl)
        self.nrUnqSamples = self.sampleSize - len(self.dupList)
        if hasattr(self, 'nrRNGs'):
            self.nrUnqSamples = self.nrUnqSamples * self.nrRNGs 
        if hasattr(self, 'propPost'):
            self.Amps = np.zeros((self.interpTime.size, self.nrUnqSamples))
            self.S    = np.zeros((self.nrUnqSamples, self.nrUnqSamples))
            self.H    = np.zeros((self.nrUnqSamples, self.nrUnqSamples))
            self.SInv = np.zeros((self.nrUnqSamples, self.nrUnqSamples))
            self.Heff = np.zeros((self.nrUnqSamples, self.nrUnqSamples))  

    def getExpectationIFG(self, metric):

        partialSum = 0.0
        partialSquareSum = 0.0

        for branch in self.branches:

            summand = branch.getMetric(metric)
            squareSummand = summand ** 2 

            for leaf in branch.leaves:
                tmpSummand = leaf.getMetric(metric) 
                summand   += tmpSummand
                squareSummand += tmpSummand ** 2 


            partialSum += summand
            partialSquareSum += squareSummand


        mTimeTrace   = partialSum / self.nrUnqSamples 
        stdTimeTrace = partialSquareSum / self.nrUnqSamples - mTimeTrace**2 
        stdTimeTrace = np.sqrt(1./(self.nrUnqSamples - 1) * stdTimeTrace)

        return mTimeTrace, stdTimeTrace
    
    def gatherGridData(self): 
        import zarr
        from src.writefiles import writeDataset
        import time
        ti = time.time() 

        PGrid = []
        Pot   = []
        CGrid = {}
        Coup  = {}
        for i in range(1,self.nrStates):
            for j in range(i+1, self.nrStates+1):
                CGrid[f"c{i}_{j}"] =  []
                Coup[f"c{i}_{j}"] =  []
                

        gridData = []
        for _branch in self.branches:
            if hasattr(self, 'parallel'):
                tmpGridData = dask.delayed(branch.fillGridData)(_branch) 
            else:
                tmpGridData = branch.fillGridData(_branch)
            gridData.append(tmpGridData)

        if hasattr(self, 'parallel'):
            gridData = dask.compute(*gridData)

        if not(hasattr(self, 'parallel')):
            for gridDatum in gridData: 
                PGrid.extend(gridDatum[0])
                Pot.extend(gridDatum[1])
                for i in range(1,self.nrStates):
                    for j in range(i+1, self.nrStates+1):
                        CGrid[f"c{i}_{j}"].extend(gridDatum[2][f"c{i}_{j}"])
                        Coup[f"c{i}_{j}"].extend(gridDatum[3][f"c{i}_{j}"])
        else:
            for gridDatum in gridData: 
                PGrid.append(gridDatum[0])
                Pot.append(gridDatum[1])
                for i in range(1,self.nrStates):
                    for j in range(i+1, self.nrStates+1):
                        CGrid[f"c{i}_{j}"].append(gridDatum[2][f"c{i}_{j}"])
                        Coup[f"c{i}_{j}"].append(gridDatum[3][f"c{i}_{j}"])

            PGrid = dask.array.concatenate(PGrid)
            Pot = dask.array.concatenate(Pot)
            for i in range(1,self.nrStates):
                for j in range(i+1, self.nrStates+1):
                    CGrid[f"c{i}_{j}"] = dask.array.concatenate(CGrid[f"c{i}_{j}"])
                    Coup[f"c{i}_{j}"] = dask.array.concatenate(Coup[f"c{i}_{j}"])

        with zarr.open('data/griddata.zarr', mode='w') as root:
            potential = root.create_group('potential')
            coupling = root.create_group('coupling')
            couplingGroups = {} 
            for i in range(1,self.nrStates):
                for j in range(i+1, self.nrStates+1):
                    couplingGroups[f"c{i}_{j}"] = \
                    coupling.create_group(f"c{i}_{j}")

            if not(hasattr(self, 'parallel')):
                PGrid = np.array(PGrid)
                Pot = np.array(Pot)

            message = writeDataset(potential, 'grid', PGrid,
                                   parallel = hasattr(self, 'parallel'))
            print(message)
            message = writeDataset(potential, 'values', Pot,
                                   parallel = hasattr(self, 'parallel'))
            print(message)

            for i in range(1,self.nrStates):
                for j in range(i+1, self.nrStates+1):
                    if not(hasattr(self, 'parallel')):
                        CGrid[f"c{i}_{j}"] = np.array(CGrid[f"c{i}_{j}"])

                    message = writeDataset(
                        couplingGroups[f"c{i}_{j}"], 'grid', 
                        CGrid[f"c{i}_{j}"], 
                        parallel = hasattr(self, 'parallel') 
                    )
                    print(message)

                    if not(hasattr(self, 'parallel')):
                        Coup[f"c{i}_{j}"] = np.array(Coup[f"c{i}_{j}"])

                    message = writeDataset(
                        couplingGroups[f"c{i}_{j}"], 'vectors', 
                        Coup[f"c{i}_{j}"],
                        parallel = hasattr(self, 'parallel') 
                    )
                    print(message)

        tf = time.time() 
        print(tf-ti)


    def propagatePost(self):
        init = True

        if hasattr(self, 'nrRNGs'):
            for branch in self.branches:
                branch.propagatePost()

        for iTime, time in enumerate(self.interpTime):
            self.buildHS(time)
            if init == True: 
                self.Amps[0,:] = self.initAmps()
            self.propagateAmps(iTime)
            
        pass
    

    def buildHS(self, time):
        self.S    = np.zeros((self.nrUnqSamples, self.nrUnqSamples))
        self.H    = np.zeros((self.nrUnqSamples, self.nrUnqSamples))
        self.SInv = np.zeros((self.nrUnqSamples, self.nrUnqSamples))
        self.Heff = np.zeros((self.nrUnqSamples, self.nrUnqSamples))  

        i = -1
        j = -1
        for branch1, branch2 in self.branchPairs:
            branch1.populateData()
            if branch1.geom != branch2.geom:
                branch2.populateData()
                j += 0
            else:
                i += 1
                j  = i

            self.S[i, j] = self.calculateS(branch1, branch2)  

        self.SInv = self.invertS()
        self.Heff = np.linalg.matmul(self.SInv, self.H) 

    def initAmps(self, time):
        raise NotImplementedError


