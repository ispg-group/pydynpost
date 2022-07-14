#!/usr/local/Cluster-Apps/python/2.7.9/bin/python
import numpy as np
import os
import math
from filesys import *
from misc import *
from parse import *
from aimsinp import *

class statePopulations(object):
    """ Class handling the calculation of adiabatic 
        electronic state populations. """
    def __init__(self, parser, cwd, dirsInCwd, psFile):
        self.prsr = parser
        self.CWD  = cwd
        self.psFile = psFile
        self.dirsInCwd = dirsInCwd

    def popInnerLoop(self, tmp_CWD, nState, nState2, nrSamples, state):
        """ Gather all populations of a given state corresponding 
            to a certain run (in the case of stochastic selection)
            and append them to the nState matrix. """                
        if self.prsr.AIMStype in ["AIMS", "ESSAIMS", "OSSAIMS","AIMSWISS"]:
            for geom in np.arange(1, self.prsr.sampleSize + 1):
                if geom in self.prsr.dupList:
                    continue
                #N_tmp = self.psFile.inputFileName(tmp_CWD, "N.dat", 
                #                                  dirType = "geom_",
                #                                  numStr = str(geom)) 
                N_tmp = self.psFile.inputFileName(tmp_CWD, "N.dat", 
                                                  dirType = self.prsr.geomDir, 
                                                  numStr = str(geom)) 
                N_tmp_dat = np.genfromtxt(N_tmp)
                time = N_tmp_dat[:,0] 
                N_tmp_state = N_tmp_dat[:,state]
                N_interp_state = np.interp(self.prsr.interpTime, time, N_tmp_state)
                nState += N_interp_state
                nState2 += N_interp_state**2
                nrSamples += 1
        elif self.prsr.AIMStype == "SINGLE":
            for sample in np.arange(1, self.prsr.sampleSize + 1):
                if sample in self.prsr.dupList:
                    continue
                N_tmp = self.psFile.inputFileName(tmp_CWD, "N.dat",
                                                  dirType = self.prsr.RNGdir, 
                                                  numStr = str(sample))
                N_tmp_dat = np.genfromtxt(N_tmp)
                time = N_tmp_dat[:,0] 
                N_tmp_state = N_tmp_dat[:,state]
                N_interp_state = np.interp(self.prsr.interpTime, time, N_tmp_state)
                nState += N_interp_state
                nState2 += N_interp_state**2
                nrSamples += 1
        else:
            dirsInTmp = getDirs(tmp_CWD)
            if dirsInTmp.size != 0:
                for geom in dirsInTmp:
                    fileName = tmp_CWD + "/" + geom + "/N.dat"
                    N_tmp_dat = np.genfromtxt(fileName)
                    time = N_tmp_dat[:,0] 
                    N_tmp_state = N_tmp_dat[:,state]
                    N_interp_state = np.interp(self.prsr.interpTime, time, N_tmp_state)
                    nState += N_interp_state
                    nState2 += N_interp_state**2
                    nrSamples += 1
            else:
                fileName = tmp_CWD + "/N.dat"
                N_tmp_dat = np.genfromtxt(fileName)
                time = N_tmp_dat[:,0] 
                N_tmp_state = N_tmp_dat[:,state]
                N_interp_state = np.interp(self.prsr.interpTime, time, N_tmp_state)
                nState += N_interp_state
                nState2 += N_interp_state**2
                nrSamples += 1

        return nState, nState2, nrSamples

    def getPopulation(self):
        """ Function that calculates the mean population of every state
            and its standard error in the full FMS Bundle. """
        for state in np.arange(1,self.prsr.nrStates+1):
            nState = np.zeros(self.prsr.interpTime.size)
            nState2 = np.zeros(self.prsr.interpTime.size)
            nrSamples = 0
            if self.dirsInCwd.size != 0:
                if self.prsr.AIMStype != "AIMS":
                    if self.prsr.AIMStype in ["ESSAIMS", "OSSAIMS","AIMSWISS"]:
                        for i in np.arange(1, self.prsr.nrRNGs + 1):
                            #tmp_CWD = self.CWD + "/rng" + str(i)
                            tmp_CWD = self.CWD + "/" + self.prsr.RNGdir + str(i)
                            nState, nState2, nrSamples = self.popInnerLoop(
                                                             tmp_CWD, nState, 
                                                             nState2, nrSamples,
                                                             state
                                                         ) 
                    elif self.prsr.AIMStype == "SINGLE":
                        nState, nState2, nrSamples = self.popInnerLoop(
                                                         self.CWD, nState,
                                                         nState2, nrSamples,
                                                         state
                                                     ) 
                    else:
                        dirsInTmp = getDirs(self.CWD)
                        for DIR in dirsInTmp:
                            tmp_CWD = self.CWD + "/" + DIR 
                            nState, nState2, nrSamples = self.popInnerLoop(
                                                             tmp_CWD, nState,
                                                             nState2, nrSamples,
                                                             state
                                                         ) 
                else:
                    tmp_CWD = self.CWD 
                    nState, nState2, nrSamples = self.popInnerLoop(
                                                    tmp_CWD, nState,
                                                    nState2, nrSamples,
                                                    state
                                                ) 
            else:
                # In case AIMS population data is only supplied as N.dat files
                for i in np.arange(1, self.prsr.sampleSize + 1):
                    if i in self.prsr.dupList:
                        continue
                    N_tmp = self.CWD + "/N_" + str(i) + ".dat" 
                    N_tmp_dat = np.genfromtxt(N_tmp)
                    N_tmp_state = N_tmp_dat[:,state]
                    time = N_tmp_dat[:,0] 
                    N_interp_state = np.interp(self.prsr.interpTime, time, N_tmp_state)
                    nState += N_interp_state
                    nState2 += N_interp_state**2

                nrSamples = (self.prsr.sampleSize - len(self.prsr.dupList))
                        

            print("")
            print("The total number of unique ICs is  \t" + str(nrSamples))

            nState_mean = nState / nrSamples
            nState_std = nState2 / nrSamples - nState_mean**2 
            nState_std = np.sqrt(1./(nrSamples - 1) * nState_std)

            print("Maximum stderr: \t" + str(np.max(nState_std)))

            stateStr = str(state)
            nState_file = "N_" + stateStr + "_" + self.prsr.AIMStype + ".dat"
            np.savetxt(nState_file, np.array([self.prsr.interpTime, nState_mean,
                                               nState_std]).T,
                       fmt="%8d %30.18e %30.18e")
