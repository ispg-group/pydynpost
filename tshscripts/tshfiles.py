#!/usr/bin/env python
import numpy as np
import os
import abc
import copy
import tempfile
import commonmethods.files as files
from abc import ABCMeta
from commonmethods.filesys import *
from commonmethods.misc import *
from commonmethods.parse import *
b2A = 0.529177249
A2b = 1./b2A
Ha2eV = 27.211386245988
#Ha2eV = 27

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

    @abc.abstractmethod
    def readKineticEnergy(self):
        pass

    @abc.abstractmethod
    def readTotalEnergy(self):
        pass

    def readPositions(self): 
        if self.prsr.code == "ABIN":
            trajFileName = 'movie.xyz'
            initFileName = 'geom'
        elif self.prsr.code == "SHARC":
            trajFileName = 'output.dat'
            initFileName = 'output.dat'

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
                    self.addPositions(currInitFileName, currTrajFileName, currFileNameRoot, nrAtoms, rngPositions)
                    geomPositions.append(rngPositions)
            else:
                currFileNameRt  = self.CWD + "/" + self.prsr.geomDir 
                currFileNameRt += str(geom)  
                currTrajFileName = currFileNameRt + "/" + trajFileName
                currInitFileName = currFileNameRt + "/" + initFileName
                self.addPositions(currInitFileName, currTrajFileName, currFileNameRt, nrAtoms, geomPositions)
            positions.append(geomPositions)
                
        return positions

class processFilesABIN(processFiles):
    def readInitState(self, fileName):
        with open(fileName + 'abin.out', 'r') as initLines:
            for initLine in initLines:
                if 'ISTATE_INIT' in initLine: 
                    initState = int(initLine.strip().split()[0][-1])

        return initState
        
    def readCurrentStates(self, pad = False): 
        fileNameRt = "pop.dat"
        currentStates = []
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
                    if pad: 
                        currState = np.interp(self.prsr.interpTime, 
                                              currTime, currState) 
                    geomCurrentStates.append(currState)
                currentStates.append(geomCurrentStates)
                     
            else:
                currFileName  = self.CWD + "/" + self.prsr.geomDir 
                currFileName += str(geom) + "/"
                restart = False
                if "restart_S0" in os.listdir(currFileName):
                    startFile = currFileName + "/restart_S0/"  
                    tZero = np.genfromtxt(startFile+"/start.dat")*0.02418884254
                    restart = True
                currFileName += fileNameRt
                currDat   = np.genfromtxt(currFileName)
                currTime  = currDat[:,0]
                currState = currDat[:,1] 
                if restart:
                    currState = currState[currTime <= tZero]
                    currState = np.append(currState, [1])
                    currTime  = currTime[currTime <= tZero]
                    currTime = np.append(currTime, [currTime[-1] + 0.01])
                else:
                    currState = currState[currTime <= self.prsr.maxTime]
                    currTime  = currTime[currTime <= self.prsr.maxTime]

                if pad: 
                    currState = np.interp(self.prsr.interpTime, 
                                          currTime, currState) 
                currentStates.append(currState)
        return currentStates
                
    def readInitStPop(self, fileName):
        initStatePop = np.zeros(np.sum(self.prsr.nrStates))
        initState = self.readInitState(fileName)
        initStatePop[initState-1] = 1.0
        return initStatePop

    def readStatePopulations(self, pad = False): 
        fileNameRt = "pop.dat"
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
                        currStatePop = np.zeros(tmpCurrStatePop.size + 1) 
                        currStatePop[0] = initStPop[state-2]   
                        currStatePop[1:] = tmpCurrStatePop  
                        currStatePop = currStatePop[currTime <= self.prsr.maxTime]
                        if pad: 
                            currStatePop = np.interp(self.prsr.interpTime, 
                                                     currTime[currTime <= self.prsr.maxTime], currStatePop) 
                        indStatePop.append(currStatePop)
                    rngStatePopulation.append(indStatePop)
                geomStatePopulation.append(rngStatePopulation)
                geomTime.append(currTime)
                     
            else:
                currFileName  = self.CWD + "/" + self.prsr.geomDir 
                currFileName += str(geom) + "/"  
                initStPop = self.readInitStPop(currFileName) 
                restart = False
                if "restart_S0" in os.listdir(currFileName):
                    startFile = currFileName + "/restart_S0/"  
                    tZero = np.genfromtxt(startFile+"/start.dat")*0.02418884254
                    restart = True
                currFileName += fileNameRt
                currDat   = np.genfromtxt(currFileName)
                currTime = np.zeros(currDat[:,0].size + 1)  
                currTime[0] = 0.0
                currTime[1:] = currDat[:,0]
                indStatePop = []
                mask = (currTime <= tZero) 
                currTime = currTime[mask] 
                currTime = np.append(currTime, [currTime[-1]+0.01]) 
                for state in np.arange(2,np.sum(self.prsr.nrStates)+2):
                    tmpCurrStatePop = currDat[:, state] 
                    currStatePop = np.zeros(currDat[:,0].size + 1)  
                    currStatePop[0] = initStPop[state-2]   
                    currStatePop[1:] = tmpCurrStatePop
                    if restart:
                        currStatePop = currStatePop[mask] 
                        if state == 2:
                            currStatePop = np.append(currStatePop, [1])
                        else:
                            currStatePop = np.append(currStatePop, [0])
                    else:
                        currStatePop = currStatePop[currTime <= self.prsr.maxTime] 
                    if pad: 
                        currStatePop = np.interp(self.prsr.interpTime, 
                                                 currTime[currTime <= self.prsr.maxTime], currStatePop) 
                    indStatePop.append(currStatePop)
                geomStatePopulation.append(indStatePop)
                geomTime.append(currTime)
                
                
        return geomStatePopulation

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
                        atomPos.extend([float(x)*b2A for x in curSpltLine[1:]])
                        curPos.append(atomPos)
                        nrLines += 1
                    else:
                        atomName = curSpltLine[0]
                        atomPos  = [atomName + str(nrLines - 1)]
                        atomPos.extend([float(x)*b2A for x in curSpltLine[1:]])
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

        

    def addPositions(self, initFileName, trajFileName, fileNameRoot, nrAtoms, outTraj):
        readTimestep = lambda a: float(a[-1]) 
        outTraj.append(self.addInitGeom(initFileName, nrAtoms))
        outTraj[-1].extend(files.addTraj(trajFileName, nrAtoms, readTimestep))
        if 'restart_S0' in os.listdir(fileNameRoot): 
            restartDir = fileNameRoot + "/restart_S0/"  
            if 'movie.xyz' in os.listdir(restartDir):
                tZero = np.genfromtxt(restartDir+"/start.dat")
                restart_S0 = files.addTraj(
                    restartDir+"/movie.xyz", 
                    nrAtoms, 
                    readTimestep,
                    tZero    
                )
                tFinal = list(dict(outTraj[-1]).keys())[-1]
                if tFinal != tZero: 
                    indTZero = list(dict(outTraj[-1]).keys()).index(float(tZero))
                    indTFinal = list(dict(restart_S0).keys()).index(tFinal)
                    outTraj[-1][indTZero+1:] = restart_S0[:indTFinal+1]
                    outTraj[-1].extend(restart_S0[indTFinal+1:])
                else:
                    outTraj[-1].extend(restart_S0)

            elif 'full' in os.listdir(restartDir):
                tZero = np.genfromtxt(restartDir+"/full/start.dat")
                restart_S0 = files.addTraj(
                    restartDir+"/full/movie.xyz", 
                    nrAtoms, 
                    readTimestep,
                    tZero
                )
                tFinal = list(dict(outTraj[-1]).keys())[-1]
                tFinalRestart = list(dict(restart_S0).keys())[-1]
                #print(tFinal, tZero, tFinalRestart)
                if ((tFinal != tZero) and (tFinal < tFinalRestart)): 
                    #print(initFileName)
                    indTZero = list(dict(outTraj[-1]).keys()).index(float(tZero))
                    indTFinal = list(dict(restart_S0).keys()).index(tFinal)
                    outTraj[-1][indTZero+1:] = restart_S0[:indTFinal+1]
                    outTraj[-1].extend(restart_S0[indTFinal+1:])
                elif ((tFinal != tZero) and (tFinal >= tFinalRestart)):
                    indTZero = list(dict(outTraj[-1]).keys()).index(float(tZero))
                    indTFinalRestart = list(dict(outTraj[-1]).keys()).index(float(tFinalRestart))
                    outTraj[-1][indTZero+1:indTFinalRestart+1] = restart_S0 
                    copyOutTraj = copy.deepcopy(outTraj[-1][:indTFinalRestart+1])
                    outTraj[-1] = copyOutTraj
                else:
                    outTraj[-1].extend(restart_S0)

    def getNrAtoms(self, fileName):
        with open(fileName, "r") as geomLines:
            return int(geomLines.readline().strip().split()[0])

    def readEnergies(self, inpDir):
        fileName = inpDir + 'PES.dat'
        potEnergies = (np.genfromtxt(fileName)[:,1:] 
                       - self.prsr.minEnergy) * Ha2eV
        fileName = inpDir + 'energies.dat'
        energyData = np.genfromtxt(fileName)
        actPotEnergy = energyData[:,1] 
        actPotEnergy = (actPotEnergy - self.prsr.minEnergy) * Ha2eV
        kinEnergy = energyData[:,2] 
        totEnergy = energyData[:,3] 
        totEnergy = (totEnergy - self.prsr.minEnergy) * Ha2eV
        return {'pot': potEnergies, 
                'kin': kinEnergy, 
                'tot': totEnergy, 
                'act': actPotEnergy,
                'time': energyData[:,0]}

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
    #def readCurrentStates(self): 
    #    pass

    #def readStatePopulations(self): 
    #    pass

    def addPositions(self, initFileName, trajFileName, fileNameRoot, nrAtoms, outTraj):
        with tempfile.NamedTemporaryFile(delete=False) as tempTrajFile:
            self.convert2XYZ(trajFileName, tempTrajFile, nrAtoms)
            tempTrajFile.close()
            readTimestep = lambda a: float(a[1]) 
            outTraj.append(files.addTraj(tempTrajFile.name, nrAtoms, readTimestep))
            #print(outTraj)


    def convert2XYZ(self, outputFile, tempTrajFile, nrAtoms):
        with open(outputFile, "r") as lines:
            readStep = False
            readGeom = False
            readElements = False
            nrLines = 0
            time = []
            traj = []
            geom = []
            elements = []
            for line in lines:
                if 'Elements' in line:
                    readElements = True
                    nrLines += 1
                    continue

                if readElements == True:
                    elements.append(line.strip()) 
                    nrLines += 1 
                    if nrLines > nrAtoms:
                        nrLines = 0 
                        readElements = False
                    continue

                if 'Step' in line:
                    readStep = True
                    continue

                if readStep == True:
                    time.append(int(line.strip())*self.prsr.step)
                    readStep = False
                    continue

                if 'Geometry' in line:
                    readGeom = True
                    nrLines += 1
                    continue

                if readGeom == True:
                    atom = [elements[nrLines-1]]
                    atom.extend([float(x)*b2A for x in line.strip().split()])
                    geom.append(atom)
                    nrLines += 1
                    if nrLines > nrAtoms:
                        nrLines = 0
                        readGeom = False
                        traj.append(geom)
                        geom = []
                    continue
                     
        for t, geom in zip(time, traj): 
            outLine = '{a:11s}'.format(a=' ')
            outLine += '{nat:3d}\n'.format(nat=nrAtoms)
            outLine += '{a:5s}t={a:7s}{t:9.5f}\n'.format(a=' ', t=t)
            for atom in geom:
                outLine += '{name:2s}'.format(name=atom[0])
                for coord in atom[1:]:
                    outLine += '{a:5s}{atm:12.10f}'.format(atm=coord, a=' ')
                outLine += '\n'
            tempTrajFile.write(str.encode(outLine))

    def getNrAtoms(self, fileName):
        nrAtoms = 0
        with open(fileName, "r") as lines:
            for line in lines: 
                if 'natom' in line:
                    nrAtoms = int(line.strip().split()[1])

        if (nrAtoms == 0):
            raise ValueError('Number of atoms is zero?!')

        return nrAtoms

    def readEnergies(self, fileName):
        fileName = fileName + 'output.dat'
        with open(fileName, "r") as lines:
            nrLines = 0
            istate = 0 
            readPotEnergies = False
            readKinEnergy = False 
            readState = False
            potEnergies = []
            kinEnergy = []
            totEnergy = []
            actPotEnergy = []

            for line in lines:
                if "Hamiltonian" in line:
                    readPotEnergies = True
                    nrLines += 1
                    potEnergy = []
                    continue

                if readPotEnergies == True:
                    nrLines += 1
                    potEnergy.append(float(line.strip().split()[istate])*Ha2eV)
                    istate += 2
                    if nrLines > np.sum(self.prsr.nrStates):
                        nrLines = 0
                        istate = 0
                        readPotEnergies = False
                        potEnergies.append(potEnergy)
                        potEnergy = []
                    continue

                if "Ekin" in line:
                    readKinEnergy = True
                    continue

                if readKinEnergy == True:
                    kinEnergy.append(float(line.strip().split()[0])*Ha2eV)
                    readKinEnergy = False
                    continue

                if " states" in line:
                    readState = True
                    continue

                if readState == True:
                    state = int(line.strip().strip()[0]) 
                    Etot = kinEnergy[-1] + potEnergies[-1][state-1]  
                    totEnergy.append(Etot)
                    actPotEnergy.append(potEnergies[-1][state-1])
                    readState = False
                    continue

        return {'pot': potEnergies, 
                'kin': kinEnergy, 
                'tot': totEnergy, 
                'act': actPotEnergy}

    def readPotEnergies(self):
        pass

    def readKineticEnergy(self):
        pass

    def readTotalEnergy(self):
        pass

    def readCurrentStates(self, pad=False): 
        fileNameRt = "output.dat"
        currentStates = []
        for geom in np.arange(1, self.prsr.sampleSize + 1):
            if geom in self.prsr.dupList:
                continue
            if self.prsr.nrRNGs != 0: 
                geomCurrentStates = []
                for rng in np.arange(1, self.prsr.nrRNGs + 1):
                    pass
            else:
                currFileName  = self.CWD + "/" + self.prsr.geomDir 
                currFileName += str(geom) + "/" + fileNameRt 
                geomCurrentStates = []
                with open(currFileName, "r") as lines:
                    readState = False
                    for line in lines:
                        if " states" in line:
                            readState = True
                            continue

                        if readState == True:
                            state = int(line.strip().strip()[0]) 
                            geomCurrentStates.append(state)
                            readState = False
                            continue

                if len(geomCurrentStates) > self.prsr.interpTime.size:
                    finalLen = len(geomCurrentStates)
                    for i in np.arange(np.abs(self.prsr.interpTime.size - finalLen)):
                       geomCurrentStates.pop()

                if pad == False:
                    currentStates.append(geomCurrentStates)
                    continue

                if len(geomCurrentStates) < self.prsr.interpTime.size:
                    finalLen = len(geomCurrentStates)
                    for i in np.arange(self.prsr.interpTime.size - finalLen):
                        geomCurrentStates.append(1)
                    currentStates.append(geomCurrentStates)
        return currentStates

    def readStatePopulations(self, pad=False): 
        fileNameRt = "output.dat"
        statePop = []
        for geom in np.arange(1, self.prsr.sampleSize + 1):
            if geom in self.prsr.dupList:
                continue
            if self.prsr.nrRNGs != 0: 
                geomCurrentStates = []
                for rng in np.arange(1, self.prsr.nrRNGs + 1):
                    pass
            else:
                currFileName  = self.CWD + "/" + self.prsr.geomDir 
                currFileName += str(geom) + "/" + fileNameRt 
                geomStatePop = []
                for i in np.arange(np.sum(self.prsr.nrStates)):
                    geomStatePop.append([])
                with open(currFileName, "r") as lines:
                    nrLines = 0
                    readStatePop = False
                    for line in lines:
                        if " Coefficients" in line:
                            readStatePop = True
                            nrLines += 1
                            continue

                        if readStatePop == True:
                            floatLine = [float(x) for x in line.strip().split()]
                            geomStatePop[nrLines-1].append(np.abs(floatLine[0]+1.j*floatLine[1])**2)
                            nrLines += 1
                            if nrLines > np.sum(self.prsr.nrStates): 
                                nrLines = 0 
                                readStatePop = False
                            continue

                if len(geomStatePop[0]) > self.prsr.interpTime.size:
                    finalLen = len(geomStatePop[0])
                    for i in np.arange(np.abs(self.prsr.interpTime.size - finalLen)):
                        for j in np.arange(np.sum(self.prsr.nrStates)):
                            geomStatePop[j].pop()

                if pad == False:
                    statePop.append(geomStatePop)
                    continue

                if len(geomStatePop[0]) < self.prsr.interpTime.size:
                    finalLen = len(geomStatePop[0])
                    for i in  np.arange(np.abs(self.prsr.interpTime.size - finalLen)):
                        geomStatePop[0].append(1)
                        for j in np.arange(1,np.sum(self.prsr.nrStates)):
                            geomStatePop[j].append(0)
                    statePop.append(geomStatePop)

        return statePop

