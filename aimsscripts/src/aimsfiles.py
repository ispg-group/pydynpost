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


    def __firstScanForce(openFile, forcesList):
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

    def __secondScanForce(openFile, forceList, tmpTSpawn, diffT):
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
        pFileType = "/forces." + str(tmpParentID) + ".xyz"
        cFileType = "/forces." + str(tmpChildID)  + ".xyz"
        parentFile = self.inputFileName(tmpCWD, pFileType) 
        childFile = self.inputFileName(tmpCWD, cFileType) 
        diffT = 0
        f = open(parentFile, "r")
        tmpForcesParent = []
        self.__firstScanForce(f, tmpForcesParent)
        # Sometimes the forces are NOT written to forces.xyz
        # so we have to look at the next nearest timestep
        if len(tmpForcesParent) == 0:
            # rewind file
            f.seek(0)
            tmpTSpawn, diffT = self.__secondScanForce(f, tmpForcesParent,
                                                      tmpTSpawn, diffT)

        # don't exactly know why this is there
        try:
            parentForces = np.genfromtxt(tmpForcesParent)[:,1:]
            parentForces = parentForces.flatten()
        except:
            parentForces = [] 

        f.close()

        f = open(childFile, "r")
        start = 0
        tmpForcesChild = []
        # only one scan necessary now, since if force is written 
        # to parents forces.xyz it is also written to that of the child. 
        self.__firstScanForce(f, tmpForcesChild)
        f.close()

        # don't exactly know why this is there
        try:
            childForces = np.genfromtxt(tmpForcesChild)[:,1:]
            childForces = childForces.flatten()
        except:
            childForces = [] 

        return parentForces, childForces, diffT

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

    def readOverlaps(self, tmpCWD, childID):
        pcOlapFile = "pcolap." + str(childID)
        overlapFile = self.inputFileName(tmpCWD, pcOlapFile)
        overlapData = np.genfromtxt(overlapFile)
        timeConnected = overlapData[:,0]
        absoluteOverlap = overlapData[:,1]**2

        return timeConnected, absoluteOverlap

    def readMomenta(self, spawnTime, ID, tmpCWD, numParticles):
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
                childID = spawn_data[3].astype(int)
                parentID = spawn_data[5].astype(int)
                nrSpawns = tSpawn.size
        else:
            #print tmpCWD 
            tSpawn   = []
            childID  = []
            parentID = [1]
            nrSpawns = 0
        
        return tSpawn, childID, parentID, nrSpawns 


    def readPositions(self, ID, tmpCWD, numParticles):
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
                        atmPos  = [atmName] 
                        for i in curline[1:]:
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

    def zeroPadArray(self, time, observable):
        # Since most TBFs are not alive for the totality
        # of the dynamics it is necessary to zero their
        # observable when their dead.  
        if time[0] > 0.0:
            newObs = np.zeros(observable.size + 20)  
            newTime = np.zeros(time.size + 20)
            newObs[20:] = observable
            newTime[20:] = time 
            newTime[:20] = np.linspace(0.0, time[0], num = 20)
        else:
            newObs = observable 
            newTime = time

        if newTime[-1] < self.prsr.maxTime:
            newerObs = np.zeros(newObs.size + 20)  
            newerTime = np.zeros(newTime.size + 20)
            newerObs[:-20] = newObs 
            newerTime[:-20] = newTime
            newerTime[-20:] = np.linspace(newTime[-1], self.prsr.maxTime, 
                                          num = 20)
        else:
            newerObs = newObs 
            newerTime = newTime  

        return newerTime, newerObs

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
                        #tmpCWD = self.CWD + "/rng" + str(rng) + "/geom_" 
                        tmpCWD  = self.CWD + "/" + self.prsr.RNGdir + str(rng) 
                        tmpCWD += "/" + self.prsr.geomDir + str(geom)
                        #tmpCWD += str(geom)
                        spawnTimes, childIDs, parentIDs, numSpawns = self.findNrSpawns(tmpCWD)
                        ampFile = "Amp.1"   
                        FGampFile = self.inputFileName(tmpCWD, ampFile)
                        FGampData = np.genfromtxt(FGampFile) 
                        if interp:
                            FGtime, FGamp = self.zeroPadArray(FGampData[:,0], 
                                                              FGampData[:,1]) 
                            interpFGamp = np.interp(self.prsr.interpTime, FGtime, 
                                                   FGamp)
                            rngTBFpop.append(interpFGamp)
                        else:
                            FGamp = FGampData[:,1]
                            rngTBFpop.append(FGamp)

                        if type(childIDs) != np.ndarray:
                            childIDs = np.array([childIDs])
                        for childID in childIDs:
                            ampFile = "Amp." + str(childID)  
                            CHampFile = self.inputFileName(tmpCWD, ampFile)
                            try:
                                CHampData = np.genfromtxt(CHampFile) 
                            except:
                                continue
                            if interp:
                                if len(list(CHampData.shape)) != 1:
                                    CHtime, CHamp = self.zeroPadArray(CHampData[:,0],
                                                                      CHampData[:,1]) 
                                    interpCHamp = np.interp(self.prsr.interpTime, CHtime, 
                                                            CHamp)
                                    rngTBFpop.append(interpCHamp)
                                else:
                                    CHtime, CHamp = self.zeroPadArray(np.array([CHampData[0]]),
                                                                      np.array([CHampData[1]])) 
                                    interpCHamp = np.interp(self.prsr.interpTime, CHtime, 
                                                            CHamp)
                                    rngTBFpop.append(interpCHamp)
                            else:
                                if len(list(CHampData.shape)) != 1:
                                    CHamp = CHampData[:,1]
                                else:
                                    CHamp = CHampData[1]
                                rngTBFpop.append(CHamp)
                                
                        geomTBFpop.append(rngTBFpop)
                else:
                    #tmpCWD = self.CWD + "/geom_" + str(geom)
                    tmpCWD = self.CWD + "/" + self.prsr.geomDir + str(geom)
                    spawnTimes, childIDs, parentIDs, numSpawns = self.findNrSpawns(tmpCWD)
                    ampFile = "Amp.1"   
                    FGampFile = self.inputFileName(tmpCWD, ampFile)
                    FGampData = np.genfromtxt(FGampFile) 
                    if interp:
                        FGtime, FGamp = self.zeroPadArray(FGampData[:,0], 
                                                          FGampData[:,1]) 
                        interpFGamp = np.interp(self.prsr.interpTime, FGtime, 
                                               FGamp)
                        geomTBFpop.append(interpFGamp)
                    else:
                        FGamp = FGampData[:,1]
                        geomTBFpop.append(FGamp)

                    if type(childIDs) != np.ndarray:
                        childIDs = np.array([childIDs])
                    for childID in childIDs:
                        ampFile = "Amp." + str(childID)  
                        CHampFile = self.inputFileName(tmpCWD, ampFile)
                        try:
                            CHampData = np.genfromtxt(CHampFile) 
                            if CHampData.shape == (4,):
                                CHampData = np.array([CHampData.tolist()]) 
                        except:
                            continue
                        if interp:
                            CHtime, CHamp = self.zeroPadArray(CHampData[:,0],
                                                              CHampData[:,1]) 
                            interpCHamp = np.interp(self.prsr.interpTime, CHtime, 
                                                    CHamp)
                            geomTBFpop.append(interpCHamp)
                        else:
                            CHamp = CHampData[:,1]
                            geomTBFpop.append(CHamp)
                        
                if len(geomTBFpop) != 0:
                    TBFpop.append(geomTBFpop)
        else:
            tmpCWD = self.CWD 
            spawnTimes, childIDs, parentIDs, numSpawns = self.findNrSpawns(tmpCWD)
            ampFile = "Amp.1"   
            FGampFile = self.inputFileName(tmpCWD, ampFile)
            FGampData = np.genfromtxt(FGampFile) 
            if interp:
                FGtime, FGamp = self.zeroPadArray(FGampData[:,0], 
                                                  FGampData[:,1]) 
                interpFGamp = np.interp(self.prsr.interpTime, FGtime, 
                                       FGamp)
                TBFpop.append(interpFGamp)
            else:
                FGamp = FGampData[:,1]
                TBFpop.append(FGamp)

            if type(childIDs) != np.ndarray:
                childIDs = np.array([childIDs])
            for childID in childIDs:
                ampFile = "Amp." + str(childID)  
                CHampFile = self.inputFileName(tmpCWD, ampFile)
                try:
                    CHampData = np.genfromtxt(CHampFile) 
                    if CHampData.shape == (4,):
                        CHampData = np.array([CHampData.tolist()]) 
                except:
                    continue
                if interp:
                    CHtime, CHamp = self.zeroPadArray(CHampData[:,0],
                                                      CHampData[:,1]) 
                    interpCHamp = np.interp(self.prsr.interpTime, CHtime, 
                                            CHamp)
                    TBFpop.append(interpCHamp)
                else:
                    CHamp = CHampData[:,1]
                    TBFpop.append(CHamp)

        return TBFpop 

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
