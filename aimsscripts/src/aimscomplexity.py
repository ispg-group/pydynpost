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

class computationalComplexity(object):
    def __init__(self, parser, cwd, dirsInCwd, psFile):
        self.prsr = parser
        self.CWD  = cwd
        self.psFile = psFile
        self.dirsInCwd = dirsInCwd

    def complexityInnerLoop(self, tmp_CWD, interpNrTBF, interpNrCalcIFG):
        # calculates the cumulative number of TBFs and theor. Nr. of ES calls
        # for every stochastic rerun
        if self.prsr.AIMStype in ["AIMS", "ESSAIMS", "OSSAIMS","AIMSWISS"]:
            for geom in np.arange(1,self.prsr.sampleSize+1):
                if geom in self.prsr.dupList:
                    continue
                if self.dirsInCwd.size != 0:
                    #fileName = self.psFile.inputFileName(tmp_CWD, "FMS.out",
                    #                                     dirType = "geom_",
                    #                                     numStr = str(geom)) 
                    fileName = self.psFile.inputFileName(tmp_CWD, "FMS.out",
                                                         dirType = self.prsr.geomDir,
                                                         numStr = str(geom)) 
                else:
                    fileName = tmp_CWD + "/FMS_" + str(geom) + ".out"
                  
                timestep, nrTBF = self.psFile.readNrTBFs(fileName)
                timestep = np.array(timestep)
                nrTBF = np.array(nrTBF)
                if nrTBF[-1] == 0:
                    nrTBF[-1] = 1
                tmpNr = np.interp(self.prsr.interpTime, timestep, nrTBF)
                interpNrCalcIFG += tmpNr*(tmpNr + 1) / 2
                interpNrTBF += tmpNr
        else:
            # Only use this option when the subdirectory structure is different
            # from the standard one
            dirsInTmp = getDirs(tmp_CWD)
            if dirsInTmp.size != 0:
                for geom in dirsInTmp:
                    fileName = tmp_CWD + "/" + geom + "/FMS.out"
                    timestep, nrTBF = self.psFile.readNrTBFs(fileName)
                    timestep = np.array(timestep)
                    nrTBF = np.array(nrTBF)
                    if nrTBF[-1] == 0:
                        nrTBF[-1] = 1
                    tmpNr = np.interp(self.prsr.interpTime, timestep, nrTBF)
                    interpNrCalcIFG += tmpNr*(tmpNr + 1) / 2
                    interpNrTBF += tmpNr
            else:
                fileName = tmp_CWD + "/FMS.out"
                timestep, nrTBF = self.psFile.readNrTBFs(fileName)
                timestep = np.array(timestep)
                nrTBF = np.array(nrTBF)
                if nrTBF[-1] == 0:
                    nrTBF[-1] = 1
                tmpNr = np.interp(self.prsr.interpTime, timestep, nrTBF)
                interpNrCalcIFG += tmpNr*(tmpNr + 1) / 2
                interpNrTBF += tmpNr

        return interpNrTBF, interpNrCalcIFG

    
    def getComplexity(self):
        # Function calculates the average number of TBFs per run and   
        # the total theor. Nr. of ES calculations needed to perform one
        # timestep of dynamics                                         
        interpNrTBF = np.zeros(self.prsr.interpTime.size) 
        interpNrCalcIFG = np.zeros(self.prsr.interpTime.size) 
        interpNrCalcTot = np.zeros(self.prsr.interpTime.size) 
        if self.prsr.AIMStype == "AIMS":
            interpNrTBF, interpCalcIFG = self.complexityInnerLoop(self.CWD,
                                                                  interpNrTBF,
                                                                  interpNrCalcIFG)
            interpNrCalcTot += interpNrTBF * (interpNrTBF + 1) / 2
            nrICs = self.prsr.sampleSize - len(self.prsr.dupList)
        elif self.prsr.AIMStype in ["ESSAIMS", "OSSAIMS", "AIMSWISS", "SWISS_dirty"]:
            tmpInterpNrTBF = 0
            for DIR in np.arange(1, self.prsr.nrRNGs + 1): 
                #tmp_CWD = self.CWD + "/rng" + str(DIR) 
                tmp_CWD = self.CWD + "/" + self.prsr.RNGdir + str(DIR) 
                tmpInterpNrTBF, interpCalcIFG = self.complexityInnerLoop(tmp_CWD,
                                                                         tmpInterpNrTBF,
                                                                         interpNrCalcIFG)
                interpNrCalcTot += tmpInterpNrTBF * (tmpInterpNrTBF + 1) / 2
                interpNrTBF += tmpInterpNrTBF
                tmpInterpNrTBF = 0
            nrICs = (self.prsr.sampleSize - len(self.prsr.dupList)) * self.prsr.nrRNGs
        else:
            tmpInterpNrTBF = 0
            nrICs = 0
            for DIR in self.dirsInCwd:
                tmp_CWD = self.CWD + "/" + DIR 
                dirsInTmp = getDirs(tmp_CWD) 
                if dirsInTmp.size == 0:
                    nrICs += 1
                else:
                    nrICs += dirsInTmp.size
                tmpInterpNrTBF, interpCalcIFG = self.complexityInnerLoop(tmp_CWD, 
                                                                         tmpInterpNrTBF,
                                                                         interpNrCalcIFG)
                interpNrCalcTot += tmpInterpNrTBF * (tmpInterpNrTBF + 1) / 2
                interpNrTBF += tmpInterpNrTBF
                tmpInterpNrTBF = 0

        interpNrTBF = interpNrTBF/nrICs
        NTBFfile = "N_TBF_" + self.prsr.AIMStype + ".dat"
        NrCalcTotFile = "N_CALLS_TOT_" + self.prsr.AIMStype + ".dat"
        NrCalcIFGFile = "N_CALLS_IFG_" + self.prsr.AIMStype + ".dat"
        saveformat = np.zeros((self.prsr.interpTime.size,2))
        saveformat[:,0] = self.prsr.interpTime 
        saveformat[:,1] = interpNrTBF 
        np.savetxt(NTBFfile, saveformat, fmt="%8d %30.18e")
        saveformat[:,1] = interpNrCalcTot
        np.savetxt(NrCalcTotFile, saveformat, fmt="%8d %30.18e")
        saveformat[:,1] = interpNrCalcIFG
        np.savetxt(NrCalcIFGFile, saveformat, fmt="%8d %30.18e")
