#!/usr/bin/env python
import numpy as np
import os
from commonmethods.filesys import *
from commonmethods.misc import *
from commonmethods.parse import *
from .aimsinp import *

class decoherenceTimes(object):
    def __init__(self, parser, cwd, psFile):
        self.prsr = parser
        self.CWD  = cwd
        self.psFile = psFile

    def innerLoopDecoherence(self, tmpCWD, ParentForces, ChildForces, alpha, tmpChildID,
                             tmpTSpawn, decoherenceTimeAIMS, decoherenceTimeModel): 
        deltaF = ParentForces - ChildForces
        deltaFSquared = deltaF ** 2
        drate = 0
        for j in np.arange(deltaFSquared.size):
            tmpAlpha = alpha[j/3] 
            drate += deltaFSquared[j]/(4 * tmpAlpha)
        timeAIMS, overlapAIMS = self.psFile.readOverlaps(tmpCWD, tmpChildID)
        overlapAIMS = overlapAIMS[timeAIMS >= tmpTSpawn]
        overlapMax = overlapAIMS[0]
        timeAIMS = timeAIMS[timeAIMS >= tmpTSpawn]
        firstTime = timeAIMS[0]
        timeAIMS = timeAIMS - tmpTSpawn
        timeInterp = np.arange(0,timeAIMS[-1],0.01)
        overlapInterp = np.interp(timeInterp, timeAIMS, overlapAIMS)
        time10 = timeInterp[overlapInterp <= np.exp(-1) * overlapMax]
        tmpDTimeAIMS = time10[0]
        for j in np.arange(1,time10.size):
            if round(time10[j] - time10[j-1]) > 1e-2:
                tmpDTimeAIMS = time10[j]
        decoherenceTimeAIMS.append(tmpDTimeAIMS)
        timeModel = np.arange(0,10*timeAIMS[-1],0.01)  
        overlapModel = np.exp(-drate * timeModel**2) * overlapMax 
        time10Model = timeModel[overlapModel <= np.exp(-1) * overlapMax] 
        modelPCFile = "ModelPC." + str(tmpChildID)
        modelSaveFile = self.psFile.inputFileName(tmpCWD, modelPCFile)
        modelFormat = np.zeros((timeModel.size/1000,2))
        modelFormat[:,0] = timeModel[0:modelFormat.shape[0]*1000:1000] + firstTime 
        modelFormat[:,1] = overlapModel[0:modelFormat.shape[0]*1000:1000]  
        np.savetxt(modelSaveFile, modelFormat, fmt="%8d %30.18e")
        decoherenceTimeModel.append(time10Model[0])
        if len(time10Model) == 0:
            print("hello, we've got a problem here")
        
        return decoherenceTimeAIMS, decoherenceTimeModel

#    def getDecoherence(self, offset, maxspawn):
    def getDecoherence(self):
        s = 0
        decoherenceTimeAIMS = []
        decoherenceTimeModel = []
        for geom in np.arange(1, self.prsr.sampleSize + 1):
            #tmpCWD = self.CWD + "/geom_" + str(geom)
            tmpCWD = self.CWD + "/" + self.prsr.geomDir + str(geom)
            if geom in self.prsr.dupList:
                continue
            try:
                tmpTSpawn, tmpChildID, tmpParentID, tmpNrSpawns = self.psFile.findNrSpawns(tmpCWD)#, maxspawn = maxspawn) 
                if any(tmpTSpawn == -10.0):
                    continue
                if tmpTSpawn.size == 0:
                    continue
            except:
                #print("geom_" + str(geom) + "did not spawn, ignoring it")
                print(self.prsr.geomDir + str(geom) + "did not spawn, ignoring it")
                continue
            #tmpTSpawn += offset
            # Find the number of particles from the input files
            controlFile = self.psFile.inputFileName(tmpCWD, "Control.dat")
            f = open(controlFile, "r")
            for line in f:
                if "NumParticles" in line:
                    numParticles = int(line.split()[0].split("=")[1])
                    break
            f.close()
            alpha = self.psFile.readWidths(tmpCWD) 
            if tmpNrSpawns > 1:
                for spawn in np.arange(tmpNrSpawns):
                    ParentForces, ChildForces, spawnCorr = self.psFile.readSpawnForces(tmpCWD, tmpTSpawn[spawn], tmpChildID[spawn], 
                                                           tmpParentID[spawn], numParticles)
                    #print ParentForces, ChildForces
                    if (len(ParentForces) == 0) or (len(ChildForces) == 0):
                        continue

                    if spawnCorr != 0:
                        tmpTSpawn[spawn] += spawnCorr

                    try:
                        decoherenceTimeAIMS, decoherenceTimeModel = self.innerLoopDecoherence(tmpCWD, ParentForces, 
                                                                             ChildForces, alpha, tmpChildID[spawn],
                                                                             tmpTSpawn[spawn], decoherenceTimeAIMS, 
                                                                             decoherenceTimeModel) 
                    except:
                        if  len(decoherenceTimeAIMS) > len(decoherenceTimeModel):
                            decoherenceTimeAIMS = decoherenceTimeAIMS[:-1]
                        elif len(decoherenceTimeAIMS) < len(decoherenceTimeModel):
                            decoherenceTimeModel = decoherenceTimeModel[:-1]
                        continue
            else:
                ParentForces, ChildForces, spawnCorr = self.psFile.readSpawnForces(tmpCWD, tmpTSpawn, tmpChildID,
                                                       tmpParentID, numParticles)
                if (len(ParentForces) == 0) or (len(ChildForces) == 0):
                    continue

                if spawnCorr != 0:
                    tmpTSpawn += spawnCorr

                try:
                    decoherenceTimeAIMS, decoherenceTimeModel = self.innerLoopDecoherence(tmpCWD, ParentForces, 
                                                                                ChildForces, alpha, tmpChildID,
                                                                                tmpTSpawn, decoherenceTimeAIMS, 
                                                                                decoherenceTimeModel) 
                except:
                    if  len(decoherenceTimeAIMS) > len(decoherenceTimeModel):
                        decoherenceTimeAIMS = decoherenceTimeAIMS[:-1]
                    elif len(decoherenceTimeAIMS) < len(decoherenceTimeModel):
                        decoherenceTimeModel = decoherenceTimeModel[:-1]
                    continue

            print(len(decoherenceTimeAIMS), len(decoherenceTimeModel))


        assert (len(decoherenceTimeAIMS) == len(decoherenceTimeModel))
        saveformat = np.zeros((len(decoherenceTimeAIMS),2))
        saveformat[:,0] = np.array(decoherenceTimeAIMS) 
        saveformat[:,1] =  np.array(decoherenceTimeModel)
#        decoherenceFile = "decoherenceTimes" + str(offset) + ".dat"
        decoherenceFile = "decoherenceTimes.dat"
        np.savetxt(decoherenceFile, saveformat, fmt="%10.5f %10.5f")
