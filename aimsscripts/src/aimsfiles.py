#!/usr/local/Cluster-Apps/python/2.7.9/bin/python
import numpy as np
from matplotlib import pyplot as plt
import os
from matplotlib import cm
import matplotlib as mpl
import math
from filesys import *
from misc import *
from parse import *
from aimsinp import *
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

    def readPositions(self, ID, tmpCWD, numParticles,
                      addAtmNames=True, bohr=False):
        posFile = "positions." + str(ID) + ".xyz"
        trajPosFile = self.inputFileName(tmpCWD, posFile) 
        merr = True 
        try:
            f = open(trajPosFile, "r")
            merr = False 
        except:
            pass

        if not(merr):
            with open(trajPosFile, "r") as posLines:
                postimes = []
                positions = []
                startRcoord = False
                lineNr = 0
                curPos = []
                curTime = 0
                for line in posLines:
                    if "Time" in line:
                        curTime = float(line.split(",")[0].split(":")[1].strip())
                        postimes.append(curTime)
                        startRcoord = True
                        lineNr += 1
                        curPos = []
                        continue

                    if (startRcoord) and (lineNr < numParticles + 1):
                        curline = line.strip().split()
                        atmName = curline[0] + str(lineNr)
                        if addAtmNames:
                            atmPos  = [atmName] 
                        else:
                            atmPos  = []
                        for i in curline[1:]:
                            if bohr:
                                atmPos.append(float(i)*A2b)
                            else:
                                atmPos.append(float(i))
                        lineNr += 1
                        curPos.append(atmPos)
                    elif (startRcoord) and (lineNr == numParticles + 1):
                        lineNr = 0
                        startRcoord = False 
                        positions.append(curPos)
                else: 
                    positions.append(curPos)

            return merr, list(zip(postimes, positions))

        return merr, [ ]

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

    def readNrTBFs(self, fileName):
        f = open(fileName, "r")
        fileContents = []
        for line in f:
            if line != "\n" and ("--Time:" in line):
                fileContents.append(line[:line.find('\n')]) 

        timestep = []  
        nrTBF = []
        for line in fileContents:
            tmp_line = line.split() 
            for i in np.arange(0,len(tmp_line)):
                if (tmp_line[i] == 'trajectories)'
                    or tmp_line[i] == 'trajectory)'):
                    nrTBF.append(int(tmp_line[i-1]))
                    break
            
            if "Centroid" not in line:
                timestep.append(float(tmp_line[1]))  

        return timestep, nrTBF

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
