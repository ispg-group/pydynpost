#!/usr/bin/env python
import numpy as np
from scipy import interpolate
import os
import time
from src.filesys import *
from src.misc import *
from src.parse import *
from src.aims.inp import *
b2A = 0.529177249
A2b = 1./b2A

class processFiles(object):
    """ 
        This class simply contains a bunch of file parsing methods for the 
        outputfiles of FMS90.
    """
    def __init__(self, parser, cwd):
        self.prsr = parser
        self.CWD  = cwd

    def inputFileName(self, tmpCWD, fileType, dirType = None, numStr = None):
        # This general function sets the absolute path of an input (rarely output) 
        # file. Since TeraChem and Molpro have different subdirectory structure
        # this is manditory. 
        fName = ""
        if self.prsr.pckg == "molpro":
            if dirType != None:
                fName  = tmpCWD + "/" + dirType + numStr + "/" + self.prsr.outputDir
                fName += "/" + fileType
            else:
                fName = tmpCWD + "/" + self.prsr.outputDir + "/" + fileType
        else:
            if dirType != None:
                fName = tmpCWD + "/" + dirType + numStr + "/" + fileType
            else:
                fName = tmpCWD + "/" + fileType
        return fName


    def __firstScanForce(self,openFile, forcesList, tmpTSpawn,
                         numParticles):
        # First scan through the parents forces.xyz file  
        start = 0
        for line in openFile:
            if (str(tmpTSpawn) in line) and ("Time" in line):
                start += 1 
                continue

            if start > 0 and start < numParticles + 1: 
                forcesList.append(line.strip())
                start += 1
            elif (start >= numParticles + 1):
                break

    def __secondScanForce(self, openFile, forceList, tmpTSpawn, diffT,
                          numParticles):
        start = 0
        for line in openFile:
            if ("Time" in line) and (start == 0):
                timestep = line.split(",")[0].split(":")[1].strip()
                diffT = float(timestep) - float(tmpTSpawn)
                if ((diffT > 0) and (diffT < 20)):
                    start += 1
                    continue
            if start > 0 and start < numParticles + 1: 
                forceList.append(line.strip())
                start += 1
            elif (start >= numParticles + 1):
                break
        
        tmpTSpawn += diffT
        return tmpTSpawn, diffT
    
    def readSpawnForces(self, tmpCWD, tmpTSpawn, tmpChildID, tmpParentID,
                        numParticles):
        pFileType = "forces." + str(tmpParentID) + ".xyz"
        cFileType = "forces." + str(tmpChildID)  + ".xyz"
        parentFile = self.inputFileName(tmpCWD, pFileType) 
        childFile = self.inputFileName(tmpCWD, cFileType) 
        diffT = 0
        f = open(parentFile, "r")
        tmpForcesParent = []
        self.__firstScanForce(f, tmpForcesParent, tmpTSpawn,
                              numParticles)
        # Sometimes the forces are NOT written to forces.xyz
        # so we have to look at the next nearest timestep
        if len(tmpForcesParent) == 0:
            # rewind file
            f.seek(0)
            tmpTSpawn, diffT = self.__secondScanForce(f, tmpForcesParent,
                                                      tmpTSpawn, diffT, 
                                                      numParticles)

        # don't exactly know why this is there
        try:
            parentForces = np.genfromtxt(tmpForcesParent)[:,1:]
            parentForces = parentForces.flatten()
        except:
            parentForces = [] 

        f.close()

        if self.prsr.pckg == "molpro":
            filesInTmpCWD = os.listdir(tmpCWD + "/" + self.prsr.outputDir)
        else:
            filesInTmpCWD = os.listdir(tmpCWD)
        if cFileType in filesInTmpCWD:
            f = open(childFile, "r")
            start = 0
            tmpForcesChild = []
            # only one scan necessary now, since if force is written 
            # to parents forces.xyz it is also written to that of the child. 
            self.__firstScanForce(f, tmpForcesChild, tmpTSpawn,
                                  numParticles)
            f.close()

            #print(tmpForcesChild)

            # don't exactly know why this is there
            try:
                childForces = np.genfromtxt(tmpForcesChild)[:,1:]
                childForces = childForces.flatten()
            except:
                childForces = [] 

            return parentForces, childForces, diffT

        else:
            return [], [], []

    def readWidths(self, tmpCWD):
        FMSFile = self.inputFileName(tmpCWD, "FMS.out")
        f = open(FMSFile, "r")
        start = 0
        tmpWidths = []
        for line in f:
            if "Particle #" in line:
                start += 1
                continue
            
            if (start > 0) and (start < 3):
                start += 1
            elif (start == 3):
                tmpWidths.append(float(line.split(":")[1].strip()))
                start = 0

        return np.array(tmpWidths) 

    def readAtomNames(self, tmpCWD):
        FMSFile = self.inputFileName(tmpCWD, "FMS.out")
        f = open(FMSFile, "r")
        start = 0
        tmpNames = []
        num = 1
        for line in f:
            if "Particle #" in line:
                start += 1
                continue
            
            if (start > 0) and (start < 4):
                start += 1
            elif (start == 4):
                tmpNames.append(line.split(":")[1].strip() + str(num))
                num += 1
                start = 0

        return np.array(tmpNames) 

    def readOverlaps(self, tmpCWD, childID):
        pcOlapFile = "pcolap." + str(childID)
        overlapFile = self.inputFileName(tmpCWD, pcOlapFile)
        overlapData = np.genfromtxt(overlapFile)
        timeConnected = overlapData[:,0]
        absoluteOverlap = overlapData[:,1]**2

        return timeConnected, absoluteOverlap

    def readSpawnMomenta(self, spawnTime, ID, tmpCWD, numParticles):
        momFile = "TrajDump." + str(ID) 
        trajDumpFile = self.inputFileName(tmpCWD, momFile) 
        trajDumpData = np.genfromtxt(trajDumpFile) 
        if spawnTime not in trajDumpData[:,0]:
            return np.array([None])
        else:
            index = np.where(trajDumpData[:,0] == spawnTime)
            index = int(index[0][0])
            return trajDumpData[index,3*numParticles+1:6*numParticles+1]

    def findNumAtoms(self, tmpCWD):
        controlFile = self.inputFileName(tmpCWD, "Control.dat")
        f = open(controlFile, "r")
        for line in f:
            if "NumParticles" in line:
                numParticles = int(line.split()[0].split("=")[1])
                break
        return numParticles

    def findNrSpawns(self, tmpCWD, maxspawn = None): 
        spawn_log = self.inputFileName(tmpCWD, "Spawn.log") 
        if self.prsr.pckg == "molpro":
            filesInTmpCWD = os.listdir(tmpCWD + "/" + self.prsr.outputDir)
        else:
            filesInTmpCWD = os.listdir(tmpCWD)
        if "Spawn.log" in filesInTmpCWD:
            spawn_data = np.genfromtxt(spawn_log)
            if spawn_data.shape != (spawn_data.size,):
                tSpawns = spawn_data[:,1]
                if maxspawn != None:
                    tSpawn = tSpawns[tSpawns <= self.prsr.maxTime]
                    childID = spawn_data[:,3][tSpawns <= self.prsr.maxTime].astype(int)
                    parentID = spawn_data[:,5][tSpawns <= self.prsr.maxTime].astype(int)
                else:
                    tSpawn = tSpawns
                    childID = spawn_data[:,3].astype(int)
                    parentID = spawn_data[:,5].astype(int)

                if tSpawn.size == 1:
                    tSpawn = tSpawn[0] 
                    childID = childID[0]
                    parentID = parentID[0]

                nrSpawns = tSpawn.size
            else:
                tSpawn = spawn_data[1]
                if maxspawn != None:
                    if tSpawn >= maxspawn:
                        tSpawn = -10.0
                childID = np.array([spawn_data[3]]).astype(int)
                parentID = np.array([spawn_data[5]]).astype(int)
                nrSpawns = tSpawn.size
        else:
            #print tmpCWD 
            tSpawn   = []
            childID  = []
            parentID = [1]
            nrSpawns = 0
        
        return tSpawn, childID, parentID, nrSpawns 

    def getTBFstate(self, tmpCWD): 
        spawn_log = self.inputFileName(tmpCWD, "Spawn.log") 
        if self.prsr.pckg == "molpro":
            filesInTmpCWD = os.listdir(tmpCWD + "/" + self.prsr.outputDir)
        else:
            filesInTmpCWD = os.listdir(tmpCWD)

        control = self.inputFileName(tmpCWD, "Control.dat") 
        with open(control, "r") as controlLines:
            for controlLine in controlLines:
                if 'InitState' in controlLine:
                    try: 
                        stateStr = controlLine.strip().split('=')[1]
                        parentSt = int(stateStr)
                    except:
                        stateStr = controlLine.strip().split('=')[1]
                        parentSt = int(stateStr.split()[0])
                    break
        
        if "Spawn.log" in filesInTmpCWD:
            spawn_data = np.genfromtxt(spawn_log)
            if spawn_data.shape != (spawn_data.size,):
                childSt = spawn_data[:,4].astype(int)
            else:
                childSt = np.array([spawn_data[3]]).astype(int)
            runningSt = np.zeros(childSt.size+1)
            runningSt[0] = parentSt
            runningSt[1:] = childSt
        else: 
            runningSt = np.array([parentSt])
        
        return runningSt.astype(int)


    def readMomenta(self, ID, tmpCWD, numParticles,
                    addAtmNames=True):
        momFile = "TrajDump." + str(ID) 
        merr = True
        if not(momFile in os.listdir(tmpCWD)):
            momFile = "momenta." + str(ID) + ".xyz" 
            try:
                f = open(momFile, "r")
                merr = False 
            except:
                pass

            if not(merr):
                with open(momFile, "r") as momLines:
                    time = []
                    momenta = []
                    startRcoord = False
                    lineNr = 0
                    curMom = []
                    curTime = 0
                    for line in momLines:
                        if "Time" in line:
                            curTime = float(line.split(",")[0].split(":")[1].strip())
                            time.append(curTime)
                            startRcoord = True
                            lineNr += 1
                            curMom = []
                            continue

                        if (startRcoord) and (lineNr < numParticles + 1):
                            curline = line.strip().split()
                            atmName = curline[0] + str(lineNr)
                            if addAtmNames: 
                                atmMom  = [atmName] 
                            else:
                                atmMom  = []
                            for i in curline[1:]:
                                atmMom.append(float(i))
                            lineNr += 1
                            curMom.append(atmMom)
                        elif (startRcoord) and (lineNr == numParticles + 1):
                            lineNr = 0
                            startRcoord = False 
                            momenta.append(curMom)
                    else: 
                        momenta.append(curMom)
        else:
            trajDumpFile = self.inputFileName(tmpCWD, momFile) 
            try: 
                trajDumpData = np.genfromtxt(trajDumpFile) 
                merr = False
            except: 
                pass

            momenta = []
            if not(merr):
                if not(len(trajDumpData.shape) == 1):
                    time = trajDumpData[:, 0]
                    for i in np.arange(time.size):
                        curMom = []
                        for j in np.arange(3*numParticles+1,6*numParticles+1,3):
                            curMom.append(trajDumpData[i,j:j+3].tolist())
                        #print len(curMom)
                        momenta.append(curMom)
                else:
                    time = trajDumpData[0]
                    time = np.array([time])
                    curMom = []
                    for j in np.arange(3*numParticles+1,6*numParticles+1,3):
                        curMom.append(trajDumpData[j:j+3].tolist())
                    momenta.append(curMom)


        if not(merr):
            #print list(zip(time, momenta)) 
            return merr, list(zip(time, momenta))
        else:
            return merr, []

    def zeroPadArray(self, time, observableRe, observableIm=[]):
        # Since most TBFs are not alive for the totality
        # of the dynamics it is necessary to zero their
        # observable when their dead.  
        if time[0] > 0.0:
            newObsRe = np.zeros(observableRe.size + 20)  
            if not(observableIm == []):
                newObsIm = np.zeros(observableRe.size + 20)  
            newTime = np.zeros(time.size + 20)
            newObsRe[20:] = observableRe
            if not(observableIm == []):
                newObsIm[20:] = observableIm
            newTime[20:] = time 
            newTime[:20] = np.linspace(0.0, time[0], num = 20)
        else:
            newObsRe = observableRe
            if not(observableIm == []):
                newObsIm = observableIm
            newTime = time

        if newTime[-1] < self.prsr.maxTime:
            newerObsRe = np.zeros(newObsRe.size + 20)  
            if not(observableIm == []):
                newerObsIm = np.zeros(newObsIm.size + 20)  
            newerTime = np.zeros(newTime.size + 20)
            newerObsRe[:-20] = newObsRe
            if not(observableIm == []):
                newerObsIm[:-20] = newObsIm
            newerTime[:-20] = newTime
            newerTime[-20:] = np.linspace(newTime[-1], self.prsr.maxTime, 
                                          num = 20)
        else:
            newerObsRe = newObsRe
            if not(observableIm == []):
                newerObsIm = newObsIm
            newerTime = newTime  

        if observableIm == []:
            return newerTime, newerObsRe
        else:
            return newerTime, newerObsRe, newerObsIm

    def addTBFpopulations(self, tmpCWD, TBFpop, ID, interp=True):
        ampFile = "Amp." + ID  
        ampFile = self.inputFileName(tmpCWD, ampFile)
        try:
            ampData = np.genfromtxt(ampFile) 
        except:
            return -1
        if interp:
            if len(list(ampData.shape)) != 1:
                if not(ID == '1'):
                    time, pop = self.zeroPadArray(ampData[:,0], 
                                                  ampData[:,1]) 
                else:
                    time = ampData[:,0]
                    pop = ampData[:,1]
                    if time[0] > 0.0:
                        time = np.zeros(pop.size + 1)
                        time[0] = 0.0
                        time[1:] = ampData[:,0]
                        pop = np.zeros(pop.size + 1)
                        pop[0] = 1.0
                        pop[1:] = ampData[:,1]
                    time, pop = self.zeroPadArray(time, 
                                                  pop)
                interpAmp = np.interp(self.prsr.interpTime, time, 
                                      pop)
            else:
                time, pop = self.zeroPadArray(np.array([ampData[0]]),
                                              np.array([ampData[1]])) 

                interpAmp = np.interp(self.prsr.interpTime, time, 
                                      pop)
            
            TBFpop.append(interpAmp)
        else:
            if len(list(ampData.shape)) != 1:
                pop = ampData[:,1]
            else:
                pop = ampData[1]
            TBFpop.append(pop)

    def getTBFpopulations(self, interp=True):
        TBFpop = []
        if hasattr(self.prsr, "sampleSize"):
            for geom in np.arange(1, self.prsr.sampleSize + 1):
                if geom in self.prsr.dupList:
                    continue
                geomTBFpop = []
                if self.prsr.AIMStype != "AIMS":
                    for rng in np.arange(1, self.prsr.nrRNGs + 1):
                        rngTBFpop = []
                        tmpCWD  = self.CWD + "/" + self.prsr.RNGdir + str(rng) 
                        tmpCWD += "/" + self.prsr.geomDir + str(geom)
                        spawnTimes, childIDs, parentIDs, numSpawns = self.findNrSpawns(tmpCWD)
                        self.addTBFpopulations(tmpCWD, rngTBFpop, str(1), interp=interp)
                        for childID in childIDs:
                            self.addTBFpopulations(tmpCWD, rngTBFpop, str(childID), interp=interp)
                        #print geom, rng, len(rngTBFpop)
                        geomTBFpop.append(rngTBFpop)
                else:
                    tmpCWD = self.CWD + "/" + self.prsr.geomDir + str(geom)
                    spawnTimes, childIDs, parentIDs, numSpawns = self.findNrSpawns(tmpCWD)
                    self.addTBFpopulations(tmpCWD, geomTBFpop, str(1), interp=interp)
                    for childID in childIDs:
                        self.addTBFpopulations(tmpCWD, geomTBFpop, str(childID), interp=interp)

                if len(geomTBFpop) != 0:
                    TBFpop.append(geomTBFpop)
        else:
            tmpCWD = self.CWD 
            spawnTimes, childIDs, parentIDs, numSpawns = self.findNrSpawns(tmpCWD)
            self.addTBFpopulations(tmpCWD, TBFpop, str(1), interp=interp)
            for childID in childIDs:
                self.addTBFpopulations(tmpCWD, TBFpop, str(childID), interp=interp)

        return TBFpop 

    def addTBFamplitude(self, tmpCWD, TBFamp, ID, interp=True):
        ampFile = "Amp." + ID  
        ampFile = self.inputFileName(tmpCWD, ampFile)
        try:
            ampData = np.genfromtxt(ampFile) 
        except:
            return -1
        if interp:
            if len(list(ampData.shape)) != 1:
                if not(ID == '1'):
                    time, ampRe, ampIm = self.zeroPadArray(ampData[:,0], 
                                                           ampData[:,2],
                                                           observableIm=
                                                           ampData[:,3]) 
                    interpAmpRe = np.interp(self.prsr.interpTime, time, 
                                            ampRe)
                    interpAmpIm = np.interp(self.prsr.interpTime, time, 
                                            ampIm)
                else:
                    time = ampData[:,0]
                    ampRe = ampData[:,2]
                    ampIm = ampData[:,3]
                    interpAmpRe = np.interp(self.prsr.interpTime, time, 
                                            ampRe)
                    interpAmpIm = np.interp(self.prsr.interpTime, time, 
                                            ampIm)
            else:
                time, ampRe, ampIm = self.zeroPadArray(np.array([ampData[0]]),
                                                       np.array([ampData[2]]),
                                                       observableIm=
                                                       np.array([ampData[3]])) 
                interpAmpRe = np.interp(self.prsr.interpTime, time, 
                                        ampRe)
                interpAmpIm = np.interp(self.prsr.interpTime, time, 
                                        ampIm)
            interpAmp = interpAmpRe + 1j * interpAmpIm
            TBFamp.append(interpAmp)
        else:
            if len(list(ampData.shape)) != 1:
                ampRe = ampData[:,2]
                ampIm = ampData[:,3]
            else:
                ampRe = ampData[2]
                ampIm = ampData[3]
            amp = ampRe + 1j * ampIm
            TBFamp.append(amp)

    def getTBFamplitude(self, interp=True):
        TBFamp = []
        if hasattr(self.prsr, "sampleSize"):
            for geom in np.arange(1, self.prsr.sampleSize + 1):
                if geom in self.prsr.dupList:
                    continue
                geomTBFamp = []
                if self.prsr.AIMStype != "AIMS":
                    for rng in np.arange(1, self.prsr.nrRNGs + 1):
                        rngTBFamp = []
                        tmpCWD  = self.CWD + "/" + self.prsr.RNGdir + str(rng) 
                        tmpCWD += "/" + self.prsr.geomDir + str(geom)
                        spawnTimes, childIDs, parentIDs, numSpawns = self.findNrSpawns(tmpCWD)
                        self.addTBFamplitude(tmpCWD, rngTBFamp, str(1), interp=interp)
                        for childID in childIDs:
                            self.addTBFamplitude(tmpCWD, rngTBFamp, str(childID), interp=interp)
                                
                        geomTBFamp.append(rngTBFamp)
                else:
                    tmpCWD = self.CWD + "/" + self.prsr.geomDir + str(geom)
                    spawnTimes, childIDs, parentIDs, numSpawns = self.findNrSpawns(tmpCWD)
                    self.addTBFamplitude(tmpCWD, geomTBFamp, str(1), interp=interp)
                    for childID in childIDs:
                        self.addTBFamplitude(tmpCWD, geomTBFamp, str(childID), interp=interp)

                if len(geomTBFamp) != 0:
                    TBFamp.append(geomTBFamp)
        else:
            tmpCWD = self.CWD 
            spawnTimes, childIDs, parentIDs, numSpawns = self.findNrSpawns(tmpCWD)
            self.addTBFamplitude(tmpCWD, TBFamp, str(1), interp=interp)
            for childID in childIDs:
                self.addTBFamplitude(tmpCWD, TBFamp, str(childID), interp=interp)

        return TBFamp 

    def getTBFphase(self, ID, tmpCWD, interp=True):
        trajDump = 'TrajDump.' + str(ID)
        trajDumpFile = self.inputFileName(tmpCWD, trajDump) 
        mErr = True
        if (trajDump in os.listdir(tmpCWD)):
            try: 
                trajDumpData = np.genfromtxt(trajDumpFile) 
                mErr = False
            except: 
                pass

            phases = []
            if not(mErr):
                if not(len(trajDumpData.shape) == 1):
                    time = trajDumpData[:, 0]
                    for i in np.arange(time.size):
                        phases.append(trajDumpData[i,-5])
                else: 
                    time = np.array([trajDumpData[0]])
                    phases.append(trajDumpData[-5])
                    

            if not(mErr):
                if (interp == True):
                    phasesInterp = np.interp(self.prsr.interpTime, time,
                                             phases) 
                    return mErr, phasesInterp
                else:
                    return mErr, list(zip(time, phases))
            else:
                return mErr, []
        else:
            phaseFile = 'Phase.' + str(ID)
            try:
                phaseData = np.genfromtxt(phaseFile)
                mErr = False
            except:
                pass
            
            if not(mErr):
                if not(len(phaseData.shape) == 1):
                    time = phaseData[:,0]
                    phases = phaseData[:,1]
                else:
                    time = np.array(phaseData[0])
                    phases = np.array(phaseData[1])

            if not(mErr):
                if (interp == True):
                    phasesInterp = np.interp(self.prsr.interpTime, time,
                                             phases) 
                    return mErr, phasesInterp
                else:
                    return mErr, list(zip(time, phases))
            else:
                return mErr, []


    def readCSThresh(self, tmpCWD):
        mErr = False
        fileName = tmpCWD + "/Control.dat" 
        with open(fileName, "r") as controlLines:
            CSThresh = 0.0
            for controlLine in controlLines: 
                if "CSThresh" in controlLine:
                    strtIndex = controlLine.index("=") + 1
                    CSThresh = float(controlLine[strtIndex:].split()[0])

        if CSThresh == 0.0:
            mErr = True

        return mErr, CSThresh

    def readCouplings(self, fileNr, tmpCWD):
        mErr = False 
        coupTime = []
        coup     = []
        try:
            fileName = tmpCWD + "/Coup." + str(fileNr)
            coupData = np.genfromtxt(fileName)
            coupTime = coupData[:, 0] 
            if self.prsr.couplingType == "coup":
                coup = coupData[:, 1 : self.prsr.nrStates + 1]
            elif self.prsr.couplingType == "coupV":
                coup = np.abs(coupData[:, self.prsr.nrStates + 1 : ])
        except:
            mErr = True
        return mErr, [coupTime, coup] 

def setCWDPath(glbl, geom, rng = None):
    if not (rng == None):
        cwdPath = glbl.CWD + "/" + glb.RNGdir + str(rng) + "/" 
        cwdPath += glbl.geomDir + str(geom)
    else:
        cwdPath = glbl.CWD + "/" + glbl.geomDir + str(geom)    

    return cwdPath

def resolveInputFilePath(glbl, cwdPath, fileName):
    if glbl.pckg == "molpro":
        filePath = cwdPath + "/" + glbl.outputDir + "/" + fileName
    else:
        filePath = cwdPath + "/" + fileName
    return filePath

def tryOpenFile(fileToOpen, fromNp = False, message = None):
    merr = False

    if fromNp == True:
        try:
            f = np.genfromtxt(fileToOpen)
        except: 
            merr = f"Error while opening {fileToOpen}."  
            if message != None:
                merr += message 

        return merr

    try:
        f = open(fileToOpen, "r")
        f.close()
    except:
        merr = f"Error while opening {fileToOpen}."  
        if message != None:
            merr += message 

    return merr
    


def readPopulation(glbl, geom, rng = None):

    cwdPath = setCWDPath(glbl, geom, rng)
    filePath = resolveInputFilePath(glbl, cwdPath, "N.dat") 
    merr = tryOpenFile(filePath, fromNp=True) 
    if merr != False:
        raise IOError(merr)
    rawData = np.genfromtxt(filePath)
    rawTime = rawData[:, 0] 
    rawPopulation = rawData[:, 1:rawData.shape[1]-1] 

    return rawTime, rawPopulation
    
def parseCheckpointFileFor(glbl, queries, geom, rng = None):
    
    searchDict = {'pgrid': ('Positions', 'vector', glbl.nrParticles), 
                  'pot': ('Energies', 'scalar', glbl.nrStates), 
                  'coup': ('coupling', 'vector', glbl.nrParticles)}

    transToQueries = {'Positions': 'pgrid', 'Energies': 'pot', 
                      'coupling': 'coup'} 

    keywords = []
    for query in queries:
        if query in searchDict.keys():
            keywords.append(searchDict[query])

    cwdPath = setCWDPath(glbl, geom, rng)
    filePath = resolveInputFilePath(glbl, cwdPath, "Checkpoint.txt") 
    outData = {'pgrid': [], 'pot': [], 'cgrid': {}, 'coup': {}}
    for i in range(1, glbl.nrStates):
        for j in range(i+1, glbl.nrStates+1):
            outData['cgrid'][f"c{i}_{j}"] = []
            outData['coup'][f"c{i}_{j}"] = []
    
    couplingStates = []

    merr = tryOpenFile(filePath) 

    with open(filePath, 'r') as checkpointFile:
        startRead = False
        numLine = 0
        currKeyword = ''
        noCouplingData = []
        gridNr = 0
        for line in checkpointFile:
            splitLine = line.strip().split()
            if len(splitLine)  == 0:
                continue
                
            for keyword in keywords:
                if not(keyword[0] in line):
                    continue

                if keyword[1] == 'vector': 
                    numLine = 0
                    numLines = keyword[2]
                    startRead = True 
                elif keyword[1] == 'scalar':
                    numLine = 0
                    numLines = 1
                    startRead = True 
                if keyword[0] == 'coupling':
                    involvedStates = [int(x) for x in splitLine[2:4]]
                tempList = []
                currKeyword = keyword[0]
                break

            if startRead == False:
                continue

            if numLine == 0:
                numLine += 1
                continue

            tempList.extend([float(x) for x in splitLine])

            numLine += 1
            if numLine <= numLines:
                continue

            if ((glbl.model == 'zero') and
                ((transToQueries[currKeyword] == 'pgrid') or 
                 (transToQueries[currKeyword] == 'coup'))): 
                tempList = [
                    tempList[i] for i in range(0,len(tempList),3)
                ]

            if (transToQueries[currKeyword] == 'coup'): 
                if np.allclose(tempList,np.zeros(len(tempList))):
                    startRead = False
                    continue

                st1 = involvedStates[0]
                st2 = involvedStates[1]
                if st1 > st2:
                    statePair = f"c{st2}_{st1}" 
                    tempList = [
                        -tempList[i] for i in range(len(tempList))
                    ]
                else:
                    statePair = f"c{st1}_{st2}" 

                outData['cgrid'][statePair].append(
                    outData['pgrid'][-1]
                )
                outData['coup'][statePair].append(tempList)
            else:
                outData[transToQueries[currKeyword]].append(tempList) 
            startRead = False

    return False, outData 

def framesGenerator_vector(glbl, trajFile, addAtmNames=True, bohr=False):
    def _framesGenerator_vector():
        with open(trajFile, "r") as trajLines:
            startRead = False
            lineNr = 0
            for line in trajLines:
                if "Time" in line:
                    currTime = float(
                        line.split(",")[0].split(":")[1].strip()
                    )
                    startRead = True
                    lineNr += 1
                    if addAtmNames:
                        currCoord = {}
                    else: 
                        currCoord = []
                    continue

                if startRead == False:
                    continue

                if lineNr == (glbl.nrParticles + 1):
                    lineNr = 0
                    startRead = False 

                    yield {
                        'time': currTime, 
                        'value': currCoord
                    }
                    continue

                currLine = line.strip().split()
                #print(currline)
                atmName = currLine[0] + str(lineNr)
                atmCoord  = []

                for i in currLine[1:]:
                    if bohr:
                        atmCoord.append(float(i)*A2b)
                    else:
                        atmCoord.append(float(i))
                lineNr += 1
                
                if addAtmNames == True:
                    currCoord[atmName] = np.array(atmCoord)
                else:
                    currCoord.append(atmCoord)

    return _framesGenerator_vector

def framesGenerator_scalar(trajFile, dtype = 'f8', dataRange = None,
                           reverse = False):
    def _framesGenerator_scalar():
        trajData = np.genfromtxt(trajFile) 
        if reverse == True:
            trajData = trajData[::-1] 
        for trajDatum in trajData: 
            if dataRange != None:  
                yield {
                    'time': trajDatum[0], 
                    'value': trajDatum[dataRange]
                }
                continue

            if dtype == 'f8':
                yield {
                    'time': trajDatum[0], 
                    'value': trajDatum[1]
                }
            elif dtype == 'c8': 
                yield {
                    'time': trajDatum[0], 
                    'value': trajDatum[2] + 1.j * trajDatum[3]
                }
    return _framesGenerator_scalar


def linInterpolate(time, prev, current, vector = True):
    if np.isclose(time, prev['time'], rtol=1e-7, atol=1e-10):
        return prev['value'] 

    if np.isclose(time, current['time'], rtol=1e-7, atol=1e-10):
        return current['value']
        
    stepSize = time - prev['time'] 
    if vector != True: 
        secant = (current['value'] - prev['value'])
        secant /= (current['time'] - prev['time']) 
        return prev[1] + secant * stepSize  

    interpVector = {} 
    for atmName in current['value'].keys():
        secant = (current['value'][atmName] - prev['value'][atmName])
        secant /= (current['time'] - prev['time']) 
        interpVector[atmName] = np.array(secant)
        
    return interpVector
        

def initialCondition(glbl, geom, rng = None, readMom = False):
    cwdPath = setCWDPath(glbl, geom, rng)
    posFile = "Geometry.dat"
    filePath = resolveInputFilePath(glbl, cwdPath, posFile) 
    merr = False 
    merr = tryOpenFile(filePath)

    if merr != False:
        return merr, [ ]

    with open(filePath, 'r') as geomFile:
        lineNr = 0
        atomNr = 0
        momNr = 0
        startRead = False
        currCoord = {}
        for line in geomFile:
            lineNr += 1
            if lineNr == 2: 
                startRead = True
                atomNr += 1
                continue
            
            if startRead == False:
                continue

            if (lineNr <= (2 + glbl.nrParticles)):
                currLine = line.strip().split()
                atmName = currLine[0] + str(atomNr)
                currCoord[atmName] = np.genfromtxt(currLine[1:]) 
                atomNr += 1

            if ((readMom == False) and 
                (lineNr == (2 + glbl.nrParticles))):
                break
            elif ((readMom == True) and 
                  (lineNr == (2 + glbl.nrParticles))): 
                momNr += 1
                continue
            elif lineNr > (3 + glbl.nrParticles):
                currLine = line.strip().split()
                for iKey, key in enumerate(currCoord.keys()):
                    if momNr == (iKey + 1):
                        currCoord[key] = np.genfromtxt(currLine)
                momNr += 1


    return merr, currCoord
   

#        if type(time) == list:
#            quantities = []
#            for _time in time:
#                while True:
#                    if not(_time in [prev['time'],current['time']]): 
#                        try:
#                            prev, current = current, next(frames) 
#                        except:
#                            break
#                    else: 
#                        quantity = linInterpolate(time, prev, current)
#                        quantities.append(quantity) 
#                        frames = _frames()
#                        break

def readTrajFile(glbl, trajFile, quantityType = 'vector', 
                 fromDump = False, dtype = 'f8', time = None,
                 reverse = False):

    if quantityType == 'vector': 
        if fromDump == False:
            frames = framesGenerator_vector(glbl, trajFile) 
    else:
        if fromDump == False:
            _frames = framesGenerator_scalar(trajFile, dtype)
            frames = _frames()
        else:
            _frames = framesGenerator_scalar(
                trajFile, dtype=dtype, fromDump=True, 
                dataRange=dataRange, reverse=reverse
            )
            frames = _frames()
    prev, current = next(frames), next(frames)  
        
    if time != None: 
        if quantityType == 'vector': 
            quantity = {
                atmName: np.zeros(3) 
                for atmName in prev['value'].keys() 
            }  
        else:
            quantity = 0 

        while True:
            if not(time in [prev['time'],current['time']]): 
                try:
                    prev, current = current, next(frames) 
                except:
                    break
            else: 
                quantity = linInterpolate(time, prev, current) 
                break

        return quantity 

    if quantityType == 'vector': 
        quantity = np.zeros(glbl.interpTime, glbl.nrParticles, 3) 
    else:
        quantity = np.zeros(glbl.interpTime)

    for iTime, time in enumerate(glbl.interpTime):
        if time > current['time']: 
            try:
                prev, current = current, next(frames) 
            except:
                break

            if quantityType == 'vector': 
                quantity[iTime, :, :] = linInterpolate(
                    time, prev, current
                )
            else:
                quantity[iTime] = linInterpolate(
                    time, prev, current
                )


    return quantity 

def readPositions(glbl, geom, rng = None, trajID = None, time = None,
                  fromDump = False, addAtmNames=True, bohr=False):

    if glbl.model == 'zero' and time == 0.0: 
        return initialCondition(glbl, geom, rng=rng)
    cwdPath = setCWDPath(glbl, geom, rng)
    posFile = "positions." + str(trajID) + ".xyz"
    filePath = resolveInputFilePath(glbl, cwdPath, posFile) 
    merr = False 
    merr = tryOpenFile(filePath)

    if merr != False:
        return merr, [ ]

    positions = readTrajFile(glbl, filePath, time = time) 
    return merr, positions

def readMomenta(glbl, geom, rng = None, trajID = None, 
                time = None, fromDump = False):

    if glbl.model == 'zero' and time == 0.0: 
        return initialCondition(glbl, geom, rng=rng, readMom=True)
    cwdPath = setCWDPath(glbl, geom, rng)
    momFile = "momenta." + str(trajID) + ".xyz"
    filePath = resolveInputFilePath(glbl, cwdPath, momFile) 
    merr = False 
    merr = tryOpenFile(filePath)

    if merr != False:
        return merr, [ ]

    momenta = readTrajFile(glbl, filePath, time = time) 
    return merr, momenta

def readPhase(glbl, geom, rng = None, trajID = None, time = None):

    if glbl.model == 'zero' and time == 0.0: 
        return False, 0.0
    cwdPath = setCWDPath(glbl, geom, rng)
    trajDumpFile = 'TrajDump.' + str(trajID)
    filePath = resolveInputFilePath(glbl, cwdPath, trajDumpFile) 
    merr = False
    if (trajDump in os.listdir(tmpCWD)):
        merr = tryOpenFile(filePath, fromNp=True)

        if merr != False:
            return merr, [ ]

        merr, phase = readTrajDumpFile(
            trajDumpData, dataRange=(-5,None), 
            quantityType='scalar', time=time 
        )
        return merr, phase

    phaseFile = 'Phase.' + str(trajID)
    filePath = resolveInputFilePath(glbl, cwdPath, phaseFile) 
    merr = tryOpenFile(
        filePath, fromNp=True, message=" Are you sure it exists?"
    )

    if merr != False:
        return merr, [ ]

    phase = readTrajFile(
        glbl, filePath, quantityType='scalar',
        dtype='f8', time=time 
    )

    return merr, phase

def readAmplitude(glbl, geom, rng = None, trajID = None, time = None,
                  reverse = False):
    def _readAmplitude(_time = None):
        if glbl.model == 'zero' and _time == 0.0: 
            if trajID == 1:
                initAmp = 1.0 + 0.j
            else:
                initAmp = 0.0 + 0.0j
            return False, initAmp 
        cwdPath = setCWDPath(glbl, geom, rng)
        ampFile = 'Amp.' + str(trajID)
        filePath = resolveInputFilePath(glbl, cwdPath, ampFile) 
        merr = False
        merr = tryOpenFile(filePath, fromNp=True)

        if merr != False:
            return merr, [ ]

        amp = readTrajFile(
            glbl, filePath, quantityType='scalar',
            dtype='c8', time=_time, reverse=reverse
        )

        return merr, amp

    if isinstance(time,np.ndarray):
        amps = []
        times = time[:]
        merr = False
        for time in times:
            merr, amp = _readAmplitude(_time=time) 
            if merr != False:
                raise IOError(merr)
            amps.append(amp)

        return merr, amps

    return _readAmplitude(_time=time) 
            

def readAmplitudeDot(glbl, geom, rng = None, trajID = None, time = None):
    cwdPath = setCWDPath(glbl, geom, rng)
    ampDotFile = 'AmpDot.' + str(trajID)
    filePath = resolveInputFilePath(glbl, cwdPath, ampDotFile) 
    merr = False
    merr = tryOpenFile(filePath, fromNp=True)

    if merr != False:
        return merr, [ ]

    ampDot = readTrajFile(
        glbl, filePath, quantityType='scalar',
        dtype='c8', time=time 
    )

    return merr, ampDot

def readStateID(glbl, geom, rng = None, trajID = None):
    cwdPath = setCWDPath(glbl, geom, rng)
    merr = False
    if trajID == 1: 
        controlFile = 'Control.dat'
        filePath = resolveInputFilePath(glbl, cwdPath, controlFile) 
        merr = tryOpenFile(filePath)
        
        if merr != False:
            return merr, None 

        with open(filePath, "r") as controlLines:
            for controlLine in controlLines:
                if not('InitState' in controlLine):
                    continue

                try: 
                    stateStr = controlLine.strip().split('=')[1]
                    stateID = int(stateStr)
                except:
                    stateStr = controlLine.strip().split('=')[1]
                    stateID = int(stateStr.split()[0])
                break

        return merr, stateID 

    spawnFile = "Spawn." + str(trajID)
    filePath = resolveInputFilePath(glbl, cwdPath, spawnFile) 
    merr = tryOpenFile(filePath)
    
    if merr != False:
        return merr, None 

    with open(filePath, "r") as spawnLines:
        for iSpawnLine, spawnLine in enumerate(spawnLines):
            if iSpawnLine == 1:
                stateID = int(spawnLine.strip().split()[3]) 
                break

    return merr, stateID

def readNrTBFs(glbl, geom, rng = None):
    cwdPath = setCWDPath(glbl, geom, rng)
    fmsFile = 'FMS.out'
    filePath = resolveInputFilePath(glbl, cwdPath, fmsFile) 
    merr = False
    merr = tryOpenFile(filePath)

    if merr != False:
        raise IOError(merr) 

    fileContents = []
    with open(filePath, 'r') as fmsLines:
        for fmsLine in fmsLines:
            if fmsLine != "\n" and ("--Time:" in fmsLine):
                fileContents.append(fmsLine[:fmsLine.find('\n')]) 

    timestep = []  
    nrTBF = []
    for fileContent in fileContents:
        tmp_line = fileContent.split() 
        for i in np.arange(0,len(tmp_line)):
            if (tmp_line[i] == 'trajectories)'
                or tmp_line[i] == 'trajectory)'):
                nrTBF.append(int(tmp_line[i-1]))
                break
        
        if "Centroid" not in fileContent:
            timestep.append(float(tmp_line[1]))  

    return timestep, nrTBF

def readWidths(glbl, geom, rng = None):
    cwdPath = setCWDPath(glbl, geom, rng)
    fmsFile = "FMS.out"
    filePath = resolveInputFilePath(glbl, cwdPath, fmsFile) 
    merr = False 
    merr = tryOpenFile(filePath)
    if merr != False:
        raise IOError(merr) 
    with open(filePath, "r") as f:
        start = 0
        tmpWidths = []
        for line in f:
            if "Particle #" in line:
                start += 1
                continue
            
            if (start > 0) and (start < 3):
                start += 1
            elif (start == 3):
                tmpWidth = float(line.split(":")[1].strip()) 
                tmpWidths.extend([tmpWidth for i in range(3)])
                start = 0

    return merr, np.array(tmpWidths) 

