#!/usr/bin/env python
import numpy as np
import os
from src.filesys import *
from src.misc import *
from src.parse import *
from .aimsinp import *

class analyzeSpawning(object):
    def __init__(self, parser, cwd, psFile):
        self.prsr = parser
        self.CWD  = cwd
        self.psFile = psFile

    def momentumSpawnError(self, maxspawn): 
        magnitudeError = []
        correlationFunction = []
        rmsd = []
        for geom in np.arange(1, self.prsr.sampleSize + 1):
            if geom in self.prsr.dupList:
                continue
            #tmpCWD = self.CWD + "/geom_" + str(geom)
            tmpCWD = self.CWD + "/" + self.prsr.geomDir + str(geom)
            spawnTimes, childIDs, parentIDs, numSpawns = self.psFile.findNrSpawns(tmpCWD, maxspawn)
            numParticles = self.psFile.findNumAtoms(tmpCWD)
            if numSpawns > 1:
                for i in np.arange(numSpawns):
                    parentMom = self.psFile.readSpawnMomenta(spawnTimes[i], parentIDs[i], tmpCWD, numParticles)
                    childMom  = self.psFile.readSpawnMomenta(spawnTimes[i], childIDs[i], tmpCWD, numParticles)
                    if (parentMom.any() == None) or (childMom.any() == None):
                        continue
                    else:
                        normParentMom = np.linalg.norm(parentMom)
                        normDiffMom   = np.linalg.norm(parentMom - childMom)
                        rmsd.append(normDiffMom/parentMom.size)
                        dotChildParent = np.dot(parentMom,childMom)
                        correlationFunction.append(dotChildParent/(np.linalg.norm(parentMom)*np.linalg.norm(childMom)))
                        magnitudeError.append(normDiffMom/normParentMom*100)
            else:
                parentMom = self.psFile.readSpawnMomenta(spawnTimes, parentIDs, tmpCWD, numParticles)
                childMom  = self.psFile.readSpawnMomenta(spawnTimes, childIDs, tmpCWD, numParticles)
                if (parentMom.any() == None) or (childMom.any() == None):
                    continue
                else:
                    normParentMom = np.linalg.norm(parentMom)
                    normDiffMom   = np.linalg.norm(parentMom - childMom)
                    dotChildParent = np.dot(parentMom,childMom)
                    rmsd.append(normDiffMom/parentMom.size)
                    correlationFunction.append(dotChildParent/(np.linalg.norm(parentMom)*np.linalg.norm(childMom)))
                    magnitudeError.append(normDiffMom/normParentMom*100)


        magnitudeError = np.array(magnitudeError)
        correlationFunction = np.array(correlationFunction)
        rmsd = np.array(rmsd)
        num = np.arange(1,magnitudeError.size+1)
        momentumMagFile = "momentumMagError.dat"
        saveformat = np.zeros((magnitudeError.size,2))
        saveformat[:,0] = num 
        saveformat[:,1] = magnitudeError 
        np.savetxt(momentumMagFile, saveformat, fmt="%8d %10.5f")
        momentumCorrFile = "momentumCorr.dat"
        saveformat = np.zeros((magnitudeError.size,2))
        saveformat[:,0] = num 
        saveformat[:,1] = correlationFunction 
        np.savetxt(momentumCorrFile, saveformat, fmt="%8d %10.5f")
        momentumRMSDFile = "momentumRMSD.dat"
        saveformat = np.zeros((magnitudeError.size,2))
        saveformat[:,0] = num 
        saveformat[:,1] = rmsd 
        np.savetxt(momentumRMSDFile, saveformat, fmt="%8d %10.5f")
