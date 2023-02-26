#!/usr/bin/env python
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
    def getNrAtoms(self, fileName):
        pass

    @abc.abstractmethod
    def readPotEnergies(self):
        pass

    def readPositions(self, trajFileName, initFileName): 
        if self.prsr.nrRNGs != 0: 
            fileName  = self.CWD + "/" + self.prsr.RNGdir
            fileName += str(1) + "/" + self.prsr.geomDir
            fileName += str(1) + "/" + initFileName 
        else:
            fileName  = self.CWD + "/" + self.prsr.geomDir
            fileName += str(1) + "/" + initFileName 
        nrAtoms = self.getNrAtoms(fileName) 
        positions = []
        for geom in np.arange(1, self.prsr.sampleSize + 1):
            #print "geom" + str(geom)
            geomPositions = []
            if geom in self.prsr.dupList:
                continue
            if self.prsr.nrRNGs != 0: 
                for rng in np.arange(1, self.prsr.nrRNGs + 1):
                    rngPositions = []
                    currFileNameRt  = self.CWD + "/" + self.prsr.RNGdir 
                    currFileNameRt += str(rng) + "/" + self.prsr.geomDir 
                    currFileNameRt += str(geom) 
                    currTrajFileName = currFileNameRt + "/" + trajFileName
                    currInitFileName = currFileNameRt + "/" + initFileName
                    self.addPositions(currInitFileName, currTrajFileName, nrAtoms, rngPositions)
                    geomPositions.append(rngPositions)
            else:
                currFileNameRt  = self.CWD + "/" + self.prsr.geomDir 
                currFileNameRt += str(geom)  
                currTrajFileName = currFileNameRt + "/" + trajFileName
                currInitFileName = currFileNameRt + "/" + initFileName
                self.addPositions(currInitFileName, currTrajFileName, nrAtoms, geomPositions)
            positions.append(geomPositions)
                
        #print len(positions[0])
        return positions

class processFilesABIN(processFiles):
    def readInitState(self, fileName):
        with open(fileName + 'abin.out', 'r') as initLines:
            for initLine in initLines:
                if 'ISTATE_INIT' in initLine: 
                    initState = int(initLine.strip().split()[-1][:-1])

        return initState
        
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
                    currFileName += str(geom) + "/" 
                    initState = self.readInitState(currFileName)
                    currFileName += fileNameRt
                    currDat   = np.genfromtxt(currFileName)
                    tmpCurrTime  = currDat[:,0]
                    tmpCurrState = currDat[:,1] 
                    currState = np.zeros(tmpCurrState.size + 1)
                    currTime = np.zeros(tmpCurrState.size + 1)
                    currState[0] = initState
                    currState[1:] = tmpCurrState
                    currTime[0] = 0.0
                    currTime[1:] = tmpCurrTime
                    currState = currState[currTime <= self.prsr.maxTime]
                    currTime = currTime[currTime <= self.prsr.maxTime]
                    geomCurrentStates.append(currState)
                currentStates.append(geomCurrentStates)
                     
            else:
                currFileName  = self.CWD + "/" + self.prsr.geomDir 
                currFileName += str(geom) + "/" + fileNameRt 
                currDat   = np.genfromtxt(currFileName)
                currTime  = currDat[:,0]
                currState = currDat[:,1] 
                currState = currState[currTime <= self.prsr.maxTime]
                currTime  = currTime[currTime <= self.prsr.maxTime]
                #interpCurrState = np.interp(self.prsr.interpTime, 
                #                            currTime, currState) 
                currentStates.append(currState)
        return currentStates
                
    def readInitStPop(self, fileName):
        initStatePop = np.zeros(self.prsr.nrStates)
        initState = self.readInitState(fileName)
        initStatePop[initState-1] = 1.0
        return initStatePop

    def readStatePopulations(self): 
        fileNameRt = "pop.dat"
        #statePopulation = []
        #time = []
        geomStatePopulation = []
        geomTime = []
        for geom in np.arange(1, self.prsr.sampleSize + 1):
            if geom in self.prsr.dupList:
                continue
            if self.prsr.nrRNGs != 0: 
                rngStatePopulation = []
                for rng in np.arange(1, self.prsr.nrRNGs + 1):
                    currFileName  = self.CWD + "/" + self.prsr.RNGdir 
                    currFileName += str(rng) + "/" + self.prsr.geomDir 
                    currFileName += str(geom) + "/" 
                    initStPop = self.readInitStPop(currFileName) 
                    currFileName += fileNameRt
                    currDat   = np.genfromtxt(currFileName)
                    tmpCurrTime  = currDat[:,0]
                    currTime = np.zeros(tmpCurrTime.size + 1)  
                    currTime[0] = 0.0
                    currTime[1:] = tmpCurrTime
                    indStatePop = []
                    for state in np.arange(2,self.prsr.nrStates+2):
                        tmpCurrStatePop = currDat[:, state] 
                        #interpCurrStatePop = np.interp(self.prsr.interpTime, 
                        #                               currTime, currStatePop) 
                        currStatePop = np.zeros(tmpCurrStatePop.size + 1) 
                        currStatePop[0] = initStPop[state-2]   
                        currStatePop[1:] = tmpCurrStatePop  
                        currStatePop = currStatePop[currTime <= self.prsr.maxTime]
                        indStatePop.append(currStatePop)
                        #print "pop", currStatePop.size
                        #if len(geomStatePopulation) > 0:
                        #    print currStatePop - geomStatePopulation[-1][state-2]  
                    #rngStatePopulation.append(currStatePop)
                    rngStatePopulation.append(indStatePop)
                geomStatePopulation.append(rngStatePopulation)
                currTime  = currTime[currTime <= self.prsr.maxTime]
                #print "time", currTime.size
                geomTime.append(currTime)
#                statePopulation.append(geomStatePopulation)
#                time.append(geomTime)
                     
            else:
                currFileName  = self.CWD + "/" + self.prsr.geomDir 
                currFileName += str(geom) + "/"  
                initStPop = self.readInitStPop(currFileName) 
                currFileName += fileNameRt
                currDat   = np.genfromtxt(currFileName)
                currTime = np.zeros(currDat[:,0].size + 1)  
                currTime[0] = 0.0
                currTime[1:] = currDat[:,0]
                indStatePop = []
                for state in np.arange(2,self.prsr.nrStates+2):
                    tmpCurrStatePop = currDat[:, state] 
                    currStatePop = np.zeros(tmpCurrStatePop.size + 1) 
                    #interpCurrStatePop = np.interp(self.prsr.interpTime, 
                    #                               currTime, currStatePop) 
                    currStatePop[0] = initStPop[state-2]   
                    currStatePop[1:] = tmpCurrStatePop
                    currStatePop = currStatePop[currTime <= self.prsr.maxTime] 
                    indStatePop.append(currStatePop)
                geomStatePopulation.append(indStatePop)
                #statePopulation.append(geomStatePopulation)
                #currTime  = currTime[currTime <= self.prsr.maxTime]
                ##print currTime.size
                currTime = currTime[currTime <= self.prsr.maxTime]
                geomTime.append(currTime)
                
                
        return geomStatePopulation, geomTime

    def addInitGeom(self, initFileName, nrAtoms):
        with open(initFileName, "r") as initGLines:
            readNrParticles = True
            startRead = False
            nrLines = 0
            posTimes = []
            positions = []
            curPos = []
            for initGLine in initGLines:
                curSpltLine = initGLine.strip().split()
                if (len(curSpltLine) == 1) and (readNrParticles == True):
                    nrParticles     = int(curSpltLine[0])
                    readNrParticles = False 
                    startRead = True
                    nrLines += 1
                    continue
                elif (len(curSpltLine) == 1): 
                    startRead = True
                    nrLines += 1
                    continue

                if (startRead and (nrLines >= 2)):
                    if (nrLines < (nrParticles + 1)):
                        atomName = curSpltLine[0]
                        atomPos  = [atomName + str(nrLines - 1)]
                        atomPos.extend([float(x)*0.529177 for x in curSpltLine[1:]])
                        curPos.append(atomPos)
                        nrLines += 1
                    else:
                        atomName = curSpltLine[0]
                        atomPos  = [atomName + str(nrLines - 1)]
                        atomPos.extend([float(x)*0.529177 for x in curSpltLine[1:]])
                        curPos.append(atomPos)
                        positions.append(curPos)
                        startRead = False
                        nrLines = 0
                        curPos = []
                else:
                    nrLines += 1
                    curTime = 0.0 
                    posTimes.append(curTime)
                    continue

        return list(zip(posTimes,positions))

        

    def addPositions(self, initFileName, trajFileName, nrAtoms, outTraj):
        readTimestep = lambda a: float(a[-1]) 
        outTraj.append(self.addInitGeom(initFileName, nrAtoms))
        #print outTraj
        outTraj[-1].extend(files.addTraj(trajFileName, nrAtoms, readTimestep))
        #print outTraj
        #print len(outTraj)

    def getNrAtoms(self, fileName):
        with open(fileName, "r") as geomLines:
            return int(geomLines.readline().strip().split()[0])

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
