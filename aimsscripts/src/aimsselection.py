#!/usr/bin/env python
import numpy as np
import os
from filesys import *
from misc import *
from parse import *
from aimsinp import *

class selectionConsistency(object):
    def __init__(self, parser, cwd, psFile):
        self.prsr = parser
        self.CWD  = cwd
        self.psFile = psFile


    def selectionInnerLoop(self, tmpCWD, successes, failures,
                           successDTime, failureDTime):
        for geom in np.arange(1, self.prsr.sampleSize + 1):
            if geom in self.prsr.dupList:
                continue
           
            #successFile = self.psFile.inputFileName(tmpCWD, "Select.log", 
            #                                 dirType = "geom_", 
            #                                 numStr = str(geom))
            #failureFile = self.psFile.inputFileName(tmpCWD, "FailSelect.log", 
            #                                 dirType = "geom_", 
            #                                 numStr = str(geom))
            successFile = self.psFile.inputFileName(tmpCWD, "Select.log", 
                                             dirType = self.prsr.geomDir, 
                                             numStr = str(geom))
            failureFile = self.psFile.inputFileName(tmpCWD, "FailSelect.log", 
                                             dirType = self.prsr.geomDir, 
                                             numStr = str(geom))
            tSpawn, childID, parentID, nrSpawns = self.psFile.findNrSpawns(
                                                      tmpCWD + '/' + self.prsr.geomDir + str(geom)
                                                  )
            #print tSpawn
            try:
                with open(successFile, "r") as successData:
                    for line in successData:
                        lineArr = line.split()
                        if "#" not in lineArr[0]:
                            selectionTime = float(lineArr[1].strip())
                            predictedOlap = float(lineArr[-2].strip())
                            actualOlap    = float(lineArr[-1].strip()) 
                            success = [selectionTime, predictedOlap,
                                       actualOlap]
                            successes.append(success)

            except:
                pass 

            try:
                with open(failureFile, "r") as successData:
                    for line in successData:
                        lineArr = line.split()
                        if "#" not in lineArr[0]:
                            selectionTime = float(lineArr[1].strip())
                            predictedOlap = float(lineArr[-2].strip())
                            actualOlap    = float(lineArr[-1].strip()) 
                            if (1 - (1.*predictedOlap/actualOlap)) > 0.2:
                                failure = [selectionTime, predictedOlap,
                                           actualOlap]
                                failures.append(failure)
            except:
                pass 



    def getSelection(self):
        assert(self.prsr.AIMStype == "AIMSWISS")
        successes = []
        failures = []
        sDTime = []
        fDTime = []
        for rng in np.arange(1, self.prsr.nrRNGs + 1):
            #tmpCWD = self.CWD + "/rng" + str(rng) 
            tmpCWD = self.CWD + "/" + self.prsr.RNGdir + str(rng) 
            self.selectionInnerLoop(tmpCWD, successes, failures,
                                    sDTime, fDTime)

        tot = len(successes) + len(failures)
        fracFail = len(failures) / (1. * tot) * 100
        mErr = False
        if fracFail >= 50:
            print("WARNING: more than 50% of selections are erroneous!")
            print("Exact percentage of erroneous selections:\t%(ffail)5.2f" %
                  {'ffail': fracFail})
            mErr = True 
        failFile = "N_F_SELECT.dat"
        saveFormat = np.zeros((len(failures), 2))
        olapRec = 0.0
        for i in np.arange(1, len(failures) + 1):
            saveFormat[i - 1, 0] = i
            olapRec += np.sqrt(failures[i - 1][1])
            saveFormat[i - 1, 1] = failures[i - 1][0]
        np.savetxt(failFile, saveFormat, fmt="%8d %30.18e")
        olapRec = olapRec / len(failures)
        if mErr:
            print("Recomended overalp threshold for OSSAIMS:\t%(orec)5.3f" %
                  {'orec': olapRec})

        successFile = "N_S_SELECT.dat"
        saveFormat = np.zeros((len(successes), 2))
        for i in np.arange(1, len(successes) + 1):
            saveFormat[i - 1, 0] = i
            saveFormat[i - 1, 1] = successes[i - 1][0]
        np.savetxt(successFile, saveFormat, fmt="%8d %30.18e")
