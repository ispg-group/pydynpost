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

    def getPopulation(self):
        """ Function that calculates the mean population of every state,
            as the fraction of trajectories on that state, and its 
            standard error via the quantum amplitudes. """
        currentStates    = self.psFile.readCurrentStates()
        statePopulations = self.psFile.readStatePopulations()
        nStateCl = np.zeros((self.prsr.interpTime.size,self.prsr.nrStates)) 
        nSamplesCl = 0
        for geom in np.arange(self.prsr.sampleSize-len(self.prsr.dupList)):
            if self.prsr.nrRNGs != 0: 
                for rng in np.arange(self.prsr.nrRNGs):
                    for currTime, currState in enumerate(currentStates[geom][rng]):
                        nStateCl[currTime, currState - 1] += 1
                    nSamplesCl += 1
            else: 
                for currTime, currState in enumerate(currentStates[geom]):
                    nStateCl[currTime, currState - 1] += 1
                nSamplesCl += 1

        nStateCl_mean = nStateCl / nSamplesCl
        for state in np.arange(1,self.prsr.nrStates+1):
            nSamples = 0
            nStateQm = [] 
            for geom in np.arange(self.prsr.sampleSize-len(self.prsr.dupList)):
                if self.prsr.nrRNGs != 0: 
                    for rng in np.arange(self.prsr.nrRNGs):
                        nStateQm.append(statePopulations[geom][rng][state-1])
                        nSamples += 1
                else: 
                    nStateQm.append(statePopulations[geom][state-1])
                    nSamples += 1

            print("")
            print("The total number of unique ICs is  \t" + str(nSamples))

            nStateQm_mean = np.zeros(self.prsr.interpTime.size)
            for i in np.arange(nSamples):
                nStateQm_mean += nStateQm[i]  
            nStateQm_mean = nStateQm_mean / nSamples

            nStateQm_stderr = np.zeros(self.prsr.interpTime.size) 
            for i in np.arange(nSamples):
                nStateQm_stderr += (nStateQm[i] - nStateQm_mean)**2

            nStateQm_stderr = np.sqrt(nStateQm_stderr/(nSamples * (nSamples-1)))
            print("Maximum stderr: \t" + str(np.max(nStateQm_stderr)))

            stateStr = str(state)
            nStateQm_file = "N_" + stateStr + "_Qm" + ".dat"
            np.savetxt(nStateQm_file, np.array([self.prsr.interpTime, nStateQm_mean,
                                                nStateQm_stderr]).T,
                       fmt="%8d %30.18e %30.18e")
            nStateCl_file = "N_" + stateStr + "_Cl" + ".dat"
            assert(nSamples == nSamplesCl)
            np.savetxt(nStateCl_file, np.array([self.prsr.interpTime, nStateCl_mean[:,state-1],
                                                nStateQm_stderr]).T,
                       fmt="%8d %30.18e %30.18e")
