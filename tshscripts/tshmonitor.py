#!/usr/bin/env python
import numpy as np
import os
import math
from commonmethods.filesys import *
from commonmethods.misc import *
from commonmethods.parse import *

class Monitor(object):
    def __init__(self, parser, cwd, psFile):
        self.prsr = parser
        self.CWD  = cwd
        self.psFile = psFile

    def getMonitor(self, observable=None, noPrint=False):
        if observable == None:
            monitor = getattr(self, self.prsr.monitorObservable)
        else:
            monitor = getattr(self, observable)
        return monitor(noPrint)


    def energies(self, noPrint):
        totEnergies = []
        varTotEnergies = []
        dE = []
        minI = self.prsr.interpTime.size - 1
        problematic = {'geom': [None], 'timestep': [None], 'deltaE': [None]}  
        for geom in np.arange(1,self.prsr.sampleSize+1):
            if geom in self.prsr.dupList:
                continue
            if self.prsr.nrRNGs != 0: 
                for rng in np.arange(1,self.prsr.nrRNGs+1):
                    fileName = self.CWD + '/' + self.prsr.RNGdir + str(rng) 
                    fileName += '/' + self.prsr.geomDir + str(geom) + '/' 
                    energyDict = self.psFile.readEnergies(fileName)
                    totEnergies.append(energyDict['tot'])
                    enPot = np.array(energyDict['pot'])
                    dE.append(np.abs(enPot[:,1] - enPot[:,0]))
                continue

            fileName = self.CWD + '/' + self.prsr.geomDir + str(geom) + '/' 
            energyDict = self.psFile.readEnergies(fileName)
            enTot = np.array(energyDict['tot'])
            enPot = np.array(energyDict['pot'])
            dE.append(np.abs(enPot[:,1] - enPot[:,0]))
            if enTot.size < (self.prsr.interpTime.size - 1):
                print(enTot.size, self.prsr.interpTime.size - 1)
                if enTot.size < minI:
                    minI = enTot.size
                problematic['geom'].append(geom-1)
                problematic['timestep'].append(enTot.size)
                problematic['deltaE'].append(dE[-1][-1])
                print('early', geom)
            totEnergies.append(enTot)
            varTotEnergies.append(enTot-enTot[0])
            #time = np.arange(0, 0.5*enTot.size, 0.5)
            if noPrint == True: 
                continue
            np.savetxt(fileName + 'totEn.dat', 
                       np.array([energyDict['time'], enTot]).T,
                       fmt="%8.2f %30.18e")
            np.savetxt(fileName + 'varTotEn.dat', 
                       np.array([energyDict['time'], enTot-enTot[0]]).T,
                       fmt="%8.2f %30.18e")

        nrSamples = (self.prsr.nrRNGs + 1) * (self.prsr.sampleSize-len(self.prsr.dupList))
        mTotEnergy = 0
        stdTotEnergy = 0
        for totEnergy in totEnergies: 
            mTotEnergy += totEnergy[:minI]
            stdTotEnergy += totEnergy[:minI]**2

        mTotEnergy = mTotEnergy / nrSamples 
        stdTotEnergy = stdTotEnergy / (nrSamples - 1) - nrSamples / (nrSamples - 1) * mTotEnergy**2
        stdTotEnergy = np.sqrt(stdTotEnergy)
        time = np.arange(0, 0.5*minI, 0.5)
        if noPrint == False: 
            np.savetxt('./totEnergy.dat', np.array([time, mTotEnergy, stdTotEnergy]).T,
                       fmt="%8.2f %30.18e %30.18e")
        maxVar = stdTotEnergy[0] 
        print(maxVar)
        for iVar, var in enumerate(varTotEnergies):
            if np.abs(var).max() > maxVar:
                for iv, v in enumerate(var):
                    if v < maxVar:
                        continue
                    break
                problematic['geom'].append(iVar)
                problematic['timestep'].append(iv)
                problematic['deltaE'].append(dE[iVar][iv])
        problematic['geom'].append(None)
        problematic['timestep'].append(None)
        problematic['deltaE'].append(None)
        print(problematic)

        return problematic

