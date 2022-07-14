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
        statePopulations, time = self.psFile.readStatePopulations()
        #print statePopulations
        nStateCl = np.zeros((len(time[0]),self.prsr.nrStates)) 
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

        nStateCl_m = nStateCl / nSamplesCl
        nStateCl_mean = np.zeros((self.prsr.interpTime.size,self.prsr.nrStates))
        for i in np.arange(self.prsr.nrStates): 
            nStateCl_mean[:,i] = np.interp(self.prsr.interpTime, 
                                           time[0], nStateCl_m[:,i]) 
        print len(statePopulations) 
        for state in np.arange(1,self.prsr.nrStates+1):
            nrSamples = 0
            #nStateQm = [] 
            nStateQm_mean = 0
            nStateQm_std = 0
            for geom in np.arange(self.prsr.sampleSize-len(self.prsr.dupList)):
                if self.prsr.nrRNGs != 0: 
                    for rng in np.arange(self.prsr.nrRNGs):
                        nStateQm_mean += statePopulations[geom][rng][state-1]
                        nStateQm_std += statePopulations[geom][rng][state-1]**2
                else: 
                    nStateQm_mean += statePopulations[geom][state-1]
                    nStateQm_std += statePopulations[geom][state-1]**2 
                    nrSamples += 1

            print("")
            print("The total number of unique ICs is  \t" + str(nrSamples))

            nStateQm_mean /= nrSamples
            nStateQm_std = nStateQm_std / nrSamples - nstateQm_mean**2
            nStateQm_std = np.sqrt(1/(nrSamples-1)*nStateQm_std)
            
            nStateQm_mean = np.interp(self.prsr.interpTime, time[0],
                                      nStateQm_mean)
            nStateQm_std = np.interp(self.prsr.interpTime, time[0],
                                     nStateQm_std)

            print("Maximum stderr: \t" + str(np.max(nStateQm_std)))

            stateStr = str(state)
            nStateQm_file = "N_" + stateStr + "_Qm" + ".dat"
            np.savetxt(nStateQm_file, np.array([self.prsr.interpTime, nStateQm_mean,
                                                nStateQm_std]).T,
                       fmt="%8.2f %30.18e %30.18e")
            nStateCl_file = "N_" + stateStr + "_Cl" + ".dat"
            assert(nSamples == nSamplesCl)
            np.savetxt(nStateCl_file, np.array([self.prsr.interpTime, nStateCl_mean[:,state-1],
                                                nStateQm_std]).T,
                       fmt="%8.2f %30.18e %30.18e")
