#!/usr/local/Cluster-Apps/python/2.7.9/bin/python
import numpy as np
import os
import abc
import files 
from abc import ABCMeta
from filesys import *
from misc import *
from parse import *

class processFiles(object):
    __metaclass__ = ABCMeta
    def __init__(self, parser, cwd):
        self.prsr = parser
        self.CWD  = cwd

    
    @abc.abstractmethod
    def readCurrentStates(self): 
        pass

    @abc.abstractmethod
    def readStatePopulations(self): 
        pass

    @abc.abstractmethod
    def addPositions(self, fileName, nrAtoms, outTraj):
        pass

    @abc.abstractmethod
    def getNrAtoms(self):
        pass

    @abc.abstractmethod
    def readPotEnergies(self):
        pass

    def readPositions(self): 
        fileNameRt = "movie.xyz"
        nrAtoms = self.getNrAtoms()
        positions = []
        for geom in np.arange(1, self.prsr.sampleSize + 1):
            #print "geom" + str(geom)
            geomPositions = []
            if geom in self.prsr.dupList:
                continue
            if self.prsr.nrRNGs != 0: 
                rngPositions = []
                for rng in np.arange(1, self.prsr.nrRNGs + 1):
                    currFileName  = self.CWD + "/" + self.prsr.RNGdir 
                    currFileName += str(rng) + "/" + self.prsr.geomDir 
                    currFileName += str(geom) + "/" + fileNameRt
                    self.addPositions(currFileName, nrAtoms, rngPositions)
                geomPositions.append(rngPositions)
            else:
                currFileName  = self.CWD + "/" + self.prsr.geomDir 
                currFileName += str(geom) + "/" + fileNameRt 
                self.addPositions(currFileName, nrAtoms, geomPositions)
            positions.append(geomPositions)
                
        return positions

class processFilesABIN(processFiles):
    def readCurrentStates(self): 
        fileNameRt = "pop.dat"
        currentStates = []
        #print  self.prsr.dupList
        for geom in np.arange(1, self.prsr.sampleSize + 1):
            if geom in self.prsr.dupList:
                continue
            if self.prsr.nrRNGs != 0: 
                geomCurrentStates = []
                for rng in np.arange(1, self.prsr.nrRNGs + 1):
                    currFileName  = self.CWD + "/" + self.prsr.RNGdir 
                    currFileName += str(rng) + "/" + self.prsr.geomDir 
                    currFileName += str(geom) + "/" + fileNameRt
                    currDat   = np.genfromtxt(currFileName)
                    currTime  = currDat[:,0]
                    currState = currDat[:,1] 
                    interpCurrState = np.interp(self.prsr.interpTime, 
                                                currTime, currState) 
                    geomCurrentStates.append(interpCurrState)
                currentStates.append(geomCurrentStates)
                     
            else:
                currFileName  = self.CWD + "/" + self.prsr.geomDir 
                currFileName += str(geom) + "/" + fileNameRt 
                currDat   = np.genfromtxt(currFileName)
                currTime  = currDat[:,0]
                currState = currDat[:,1] 
                interpCurrState = np.interp(self.prsr.interpTime, 
                                            currTime, currState) 
                currentStates.append(interpCurrState)
        return currentStates
                

    def readStatePopulations(self): 
        fileNameRt = "pop.dat"
        statePopulation = []
        for geom in np.arange(1, self.prsr.sampleSize + 1):
            geomStatePopulation = []
            if geom in self.prsr.dupList:
                continue
            if self.prsr.nrRNGs != 0: 
                rngStatePopulation = []
                for rng in np.arange(1, self.prsr.nrRNGs + 1):
                    currFileName  = self.CWD + "/" + self.prsr.RNGdir 
                    currFileName += str(rng) + "/" + self.prsr.geomDir 
                    currFileName += str(geom) + "/" + fileNameRt
                    currDat   = np.genfromtxt(currFileName)
                    currTime  = currDat[:,0]
                    for state in np.arange(2,self.prsr.nrStates+2):
                        currStatePop = currDat[:, state] 
                        interpCurrStatePop = np.interp(self.prsr.interpTime, 
                                                       currTime, currStatePop) 
                        rngStatePopulation.append(interpCurrStatePop)
                    geomStatePopulation.append(rngStatePopulation)
                statePopulation.append(geomStatePopulation)
                     
            else:
                currFileName  = self.CWD + "/" + self.prsr.geomDir 
                currFileName += str(geom) + "/" + fileNameRt 
                currDat   = np.genfromtxt(currFileName)
                currTime  = currDat[:,0]
                for state in np.arange(2,self.prsr.nrStates+2):
                    currStatePop = currDat[:, state] 
                    interpCurrStatePop = np.interp(self.prsr.interpTime, 
                                                   currTime, currStatePop) 
                    geomStatePopulation.append(interpCurrStatePop)
                statePopulation.append(geomStatePopulation)
                
                
        return statePopulation

    def addPositions(self, fileName, nrAtoms, outTraj):
        readTimestep = lambda a: float(a[-1]) 
        outTraj.append(files.addTraj(fileName, nrAtoms, readTimestep, outTraj))
        print len(outTraj)

    def getNrAtoms(self):
        pass

    def readPotEnergies(self):
        pass

class processFilesNEWTONX(processFiles):
    def readCurrentStates(self): 
        pass

    def readStatePopulations(self): 
        pass


    def readPositions(self): 
        pass

    def readPotEnergies(self):
        pass

class processFilesSHARC(processFiles):
    def readCurrentStates(self): 
        pass

    def readStatePopulations(self): 
        pass

    def readPositions(self): 
        pass

    def readPotEnergies(self):
        pass
