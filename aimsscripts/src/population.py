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

class statePopulations(object):
    """ Class handling the calculation of adiabatic 
        electronic state populations. """
    def __init__(self, parser, cwd, dirsInCwd, psFile):
        self.prsr = parser
        self.CWD  = cwd
        self.psFile = psFile
        self.dirsInCwd = dirsInCwd

    def popInnerLoop(self, tmp_CWD, N_state, state):
        """ Gather all populations of a given state corresponding 
            to a certain run (in the case of stochastic selection)
            and append them to the N_state matrix. """                
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
                N_state.append(N_interp_state)
        elif self.prsr.AIMStype == "SINGLE":
            for sample in np.arange(1, self.prsr.sampleSize + 1):
                if sample in self.prsr.dupList:
                    continue
                #N_tmp = self.psFile.inputFileName(tmp_CWD, "N.dat",
                #                                  dirType = "rng", 
                #                                  numStr = str(sample))
                N_tmp = self.psFile.inputFileName(tmp_CWD, "N.dat",
                                                  dirType = self.prsr.RNGdir, 
                                                  numStr = str(sample))
                N_tmp_dat = np.genfromtxt(N_tmp)
                time = N_tmp_dat[:,0] 
                N_tmp_state = N_tmp_dat[:,state]
                N_interp_state = np.interp(self.prsr.interpTime, time, N_tmp_state)
                N_state.append(N_interp_state)
        else:
            dirsInTmp = getDirs(tmp_CWD)
            if dirsInTmp.size != 0:
                for geom in dirsInTmp:
                    fileName = tmp_CWD + "/" + geom + "/N.dat"
                    N_tmp_dat = np.genfromtxt(fileName)
                    time = N_tmp_dat[:,0] 
                    N_tmp_state = N_tmp_dat[:,state]
                    N_interp_state = np.interp(self.prsr.interpTime, time, N_tmp_state)
                    N_state.append(N_interp_state)
            else:
                fileName = tmp_CWD + "/N.dat"
                N_tmp_dat = np.genfromtxt(fileName)
                time = N_tmp_dat[:,0] 
                N_tmp_state = N_tmp_dat[:,state]
                N_interp_state = np.interp(self.prsr.interpTime, time, N_tmp_state)
                N_state.append(N_interp_state)

        return N_state 

    def getPopulation(self):
        """ Function that calculates the mean population of every state
            and its standard error in the full FMS Bundle. """
        for state in np.arange(1,self.prsr.nrStates+1):
            N_state = []
            if self.dirsInCwd.size != 0:
                if self.prsr.AIMStype != "AIMS":
                    if self.prsr.AIMStype in ["ESSAIMS", "OSSAIMS","AIMSWISS"]:
                        for i in np.arange(1, self.prsr.nrRNGs + 1):
                            #tmp_CWD = self.CWD + "/rng" + str(i)
                            tmp_CWD = self.CWD + "/" + self.prsr.RNGdir + str(i)
                            N_state = self.popInnerLoop(tmp_CWD, N_state, state) 
                        N_samples = len(N_state)
                    elif self.prsr.AIMStype == "SINGLE":
                        N_state = self.popInnerLoop(self.CWD, N_state, state) 
                        N_samples = len(N_state)
                    else:
                        dirsInTmp = getDirs(self.CWD)
                        for DIR in dirsInTmp:
                            tmp_CWD = self.CWD + "/" + DIR 
                            N_state = self.popInnerLoop(tmp_CWD, N_state, state) 
                        N_samples = len(N_state)
                else:
                    tmp_CWD = self.CWD 
                    N_state = self.popInnerLoop(tmp_CWD, N_state, state) 
                    N_samples = (self.prsr.sampleSize - len(self.prsr.dupList))
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
                    N_state.append(N_interp_state)

                N_samples = (self.prsr.sampleSize - len(self.prsr.dupList))
                        

            print("")
            print("The total number of unique ICs is  \t" + str(N_samples))

            N_state_cum = np.zeros((self.prsr.interpTime.size,)) 
            for i in np.arange(N_samples):
                N_state_cum += N_state[i]
                
            N_state_mean =  N_state_cum/N_samples

            N_state_stderr = np.zeros((self.prsr.interpTime.size,)) 
            for i in np.arange(N_samples):
                N_state_stderr += (N_state[i] - N_state_mean)**2

            N_state_stderr = np.sqrt(N_state_stderr/(N_samples * (N_samples-1)))
            print("Maximum stderr: \t" + str(np.max(N_state_stderr)))

            stateStr = str(state)
            N_state_file = "N_" + stateStr + "_" + self.prsr.AIMStype + ".dat"
            np.savetxt(N_state_file, np.array([self.prsr.interpTime, N_state_mean,
                                               N_state_stderr]).T,
                       fmt="%8d %30.18e %30.18e")
