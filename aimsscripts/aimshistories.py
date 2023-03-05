#!/usr/bin/env python
import numpy as np
import os
from commonmethods.filesys import *
from commonmethods.misc import *
from commonmethods.parse import *
from .aimsinp import *

class historyStatistics(object):
    def __init__(self, parser, cwd, psFile):
        self.prsr = parser
        self.CWD  = cwd
        self.psFile = psFile

    def histInnerLoop(self, geom):
        historyStrings = []
        historyProbs   = []
        for rng in np.arange(1, self.prsr.nrRNGs + 1):
            historyString = ""
            #tmpCWD = self.prsr.CWD + "/rng" + str(rng) + "/" + geom
            tmpCWD = self.prsr.CWD + "/" + self.prsr.RNGdir 
            tmpCWD += str(rng) + "/" + geom
            FMSFile = self.psFile.inputFileName(tmpCWD, "FMS.out")
            with open(FMSFile, "r") as FMSData:
                startr = False
                starts = False
                historyProb = 1
                countl = 0
                counts = 0
                for line in FMSData:
                    if "Time to select!" in line:
                        startr = True
                    if "Block" in line:
                        starts = True
                    if startr:
                        countl += 1
                    if starts:
                        counts += 1
                        if "*" in line:
                            historyProb *= float(line.strip().split()[2])
                            counts = 0
                            starts = False
                    if (startr) and (countl == 5):
                        line1 = line.strip().split()
                        num = [x for x in line1 if isint(x)]
                        assert(len(num) == 2)
                        child  = num[0] 
                        parent = num[1] 
                    if (startr) and ("Force Killing" in line):
                        forceKillLine = line.strip().split()
                        killedTraj = forceKillLine[4] 
                        if killedTraj == child:
                            countl = 0 
                            startr = False
                            historyString += "1"
                        elif killedTraj == parent:
                            countl = 0 
                            startr = False
                            historyString += "0"

            historyStrings.append(historyString)
            historyProbs.append(historyProb)

        print(historyStrings)
        print(dict(zip(historyStrings,historyProbs)))
        print(abs(1 - sum(list(set(historyProbs)))))
        return abs(1 - sum(list(set(historyProbs))))


                    

    def getHistories(self):
        assert(self.prsr.AIMStype == "SWISS")
        consistencyErr = 0
        for geom in np.arange(1, self.prsr.sampleSize + 1):
            if geom in self.prsr.dupList:
                continue
            #geomName = "geom_" + str(geom)
            geomName = self.prsr.geomDir + str(geom)
            consistencyErr += self.histInnerLoop(geomName)

        print(consistencyErr/self.prsr.sampleSize)
