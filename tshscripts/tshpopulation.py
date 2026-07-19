#!/usr/bin/env python
import numpy as np
import os
import math
from commonmethods.filesys import *
from commonmethods.misc import *
from commonmethods.parse import *
import tshscripts.tshmonitor as monitor

class statePopulations(object):
    """ Class handling the calculation of adiabatic 
        electronic state populations. """
    def __init__(self, parser, cwd, dirsInCwd, psFile):
        self.prsr = parser
        self.CWD  = cwd
        self.psFile = psFile
        self.dirsInCwd = dirsInCwd
        self.monitor = monitor.Monitor(parser, cwd, psFile)

    def getPopulation(self):
        """ Function that calculates the mean population of every state,
            as the fraction of trajectories on that state, and its 
            standard error via the quantum amplitudes. """
        problematicCases = self.monitor.getMonitor('energies', noPrint=True)
        print(len(problematicCases['timestep']))
        currentStates = self.psFile.readCurrentStates(pad=True)
        statePopulations = self.psFile.readStatePopulations(pad=True)
        nStateCl = np.zeros((self.prsr.interpTime.size,np.sum(self.prsr.nrStates))) 
        nSamplesCl = 0
        pdaw = False
        if 'pdaw.dat' in os.listdir(self.CWD):  
            pdaw = True 
            redWeights = np.genfromtxt('pdaw.dat')[:,1] 
            weights = []
            for i in range(1,self.prsr.sampleSize+1): 
                if i in  self.prsr.dupList:
                    print(i)
                    continue
                weights.append(redWeights[i-1]) 
            weights = np.array(weights)
            weights = weights/np.sum(weights)
            nState_pdaw = np.zeros((self.prsr.interpTime.size,np.sum(self.prsr.nrStates)))  

        for geom in np.arange(self.prsr.sampleSize-len(self.prsr.dupList)):
            if self.prsr.nrRNGs != 0: 
                for rng in np.arange(self.prsr.nrRNGs):
                    for currTime, currState in enumerate(currentStates[geom][rng]):
                        nStateCl[currTime, int(currState) - 1] += 1
                    nSamplesCl += 1
            else: 
                if (geom in problematicCases['geom']): 
                    nrProb = problematicCases['geom'].index(geom) 
                else: 
                    nrProb = None 

                if (nrProb == None):
                    pass
                elif  problematicCases['timestep'][nrProb] > currentStates[geom].size:
                    print(0, geom+1, problematicCases['geom'][nrProb]+1)
                    nrProb += 1 
                #elif ((geom in problematicCases['geom'])
                #        and (currentStates[geom][problematicCases['timestep'][nrProb]] != 1)): 
                #    print(1, geom+1, problematicCases['geom'][nrProb]+1)
                #    #print(currentStates[geom][problematicCases['timestep'][nrProb]])
                #    nrProb += 1 
                #    continue
                #elif ((geom in problematicCases['geom']) 
                #        and (currentStates[geom][problematicCases['timestep'][nrProb]] == 1) 
                #        and (problematicCases['deltaE'][nrProb] < self.prsr.gapThresh)):
                #    nrProb += 1 
                #    print(2, geom+1, problematicCases['geom'][nrProb]+1)
                #    continue
                else:
                    nrProb += 1

                #print(weights[geom])
                for currTime, currState in enumerate(currentStates[geom]):
                    nStateCl[currTime, int(currState) - 1] += 1
                    if pdaw == True:
                        nState_pdaw[currTime, int(currState) - 1] += weights[geom]
                nSamplesCl += 1

        nStateCl_m = nStateCl / nSamplesCl
        nState_pdaw = nState_pdaw / np.sum(weights)
        
        #nStateCl_mean = np.zeros((self.prsr.interpTime.size,np.sum(self.prsr.nrStates)))
        #for i in np.arange(np.sum(self.prsr.nrStates)): 
        #    nStateCl_mean[:,i] = np.interp(self.prsr.interpTime, 
        #                                   time[0], nStateCl_m[:,i]) 
        for state in np.arange(1,np.sum(self.prsr.nrStates)+1):
            nrProb = 1
            nSamples = 0
            nStateQm = [] 
            if pdaw: 
                nState_Qm_pdaw = []  
                nState_Qm_pdaw_sq = []  
            for geom in np.arange(self.prsr.sampleSize-len(self.prsr.dupList)):
                if self.prsr.nrRNGs != 0: 
                    for rng in np.arange(self.prsr.nrRNGs):
                        nStateQm.append(statePopulations[geom][rng][state-1])
                        nSamples += 1
                else: 
                    #print(problematicCases['timestep'][nrProb],currentStates[geom].size)
                    if (geom in problematicCases['geom']): 
                        nrProb = problematicCases['geom'].index(geom) 
                    else: 
                        nrProb = None 
                    if (nrProb == None):
                        pass
                    elif  problematicCases['timestep'][nrProb] > currentStates[geom].size:
                        nrProb += 1 
                        pass
                    #elif ((geom in problematicCases['geom'])
                    #        and (currentStates[geom][problematicCases['timestep'][nrProb]] != 1)):
                    #    nrProb += 1  
                    #    continue
                    #elif ((geom in problematicCases['geom']) 
                    #        and (currentStates[geom][problematicCases['timestep'][nrProb]] == 1) 
                    #        and (problematicCases['deltaE'][nrProb] < self.prsr.gapThresh)):
                    #    nrProb += 1  
                    #    continue
                    else:
                        nrProb += 1
                        
                    nStateQm.append(statePopulations[geom][state-1])
                    if pdaw:
                        nState_Qm_pdaw.append(statePopulations[geom][state-1])
                    nSamples += 1

            print("")
            print("The total number of unique ICs is  \t" + str(nSamples))

            nStateQm_m = np.zeros(nStateCl.shape[0])
            if pdaw: 
                nState_Qm_pdaw_m = np.zeros(nStateCl.shape[0])
            for i in np.arange(nSamples):
                nStateQm_m += nStateQm[i]  
                if pdaw: 
                    nState_Qm_pdaw_m += nState_Qm_pdaw[i] * weights[i] 
            nStateQm_m = nStateQm_m / nSamples
            if pdaw:
                nState_Qm_pdaw_m = nState_Qm_pdaw_m / np.sum(weights)
            #nStateQm_mean = np.interp(self.prsr.interpTime, time[0],
            #                          nStateQm_m)

            nStateQm_std = np.zeros(nStateCl.shape[0]) 
            if pdaw:
                nState_Qm_pdaw_std = np.zeros(nStateCl.shape[0])
            for i in np.arange(nSamples):
                nStateQm_std += (nStateQm[i] - nStateQm_m)**2
                if pdaw:
                    nState_Qm_pdaw_std += (nState_Qm_pdaw[i] - nState_Qm_pdaw_m)**2 * weights[i]

            nStateQm_std = np.sqrt(nStateQm_std/(nSamples * (nSamples-1)))
            if pdaw:
                nState_Qm_pdaw_std = np.sqrt(nState_Qm_pdaw_std/nSamples)
            #nStateQm_stderr = np.interp(self.prsr.interpTime, time[0],
            #                            nStateQm_std)
            print("Maximum stderr: \t" + str(np.max(nStateQm_std)))

            stateStr = str(state)
            nStateQm_file = "N_" + stateStr + "_Qm" + ".dat"
            np.savetxt(nStateQm_file, np.array([self.prsr.interpTime, nStateQm_m,
                                                nStateQm_std]).T,
                       fmt="%8.2f %30.18e %30.18e")
            nStateCl_file = "N_" + stateStr + "_Cl" + ".dat"
            assert(nSamples == nSamplesCl)
            np.savetxt(nStateCl_file, np.array([self.prsr.interpTime, nStateCl_m[:,state-1],
                                                nStateQm_std]).T,
                       fmt="%8.2f %30.18e %30.18e")
            nStateCl_file = "N_" + stateStr + "_PDAW" + ".dat"
            assert(nSamples == nSamplesCl)
            np.savetxt(nStateCl_file, np.array([self.prsr.interpTime, nState_pdaw[:,state-1],
                                                nState_Qm_pdaw_std]).T,
                       fmt="%8.2f %30.18e %30.18e")
            
