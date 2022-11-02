#!/usr/bin/env python
import numpy as np
from matplotlib import pyplot as plt
import abc
from abc import ABCMeta
import os
from matplotlib import cm
import matplotlib as mpl
import math
import copy
from filesys import *
from misc import *
from parse import *
from aimsinp import *
from aimssave import *
from internals import *
from aimsexpectationvalues import *

class molpop(object):
    """
        This class handles the calculation of internals  
        when there are many equivalent bonds that can 
        break. This is done in the same way as described 
        in: 
            Crespo-Otero, R. and Barbatti, M.; J. Chem. Phys. 
            134, 164305 (2011)
    """
    def __init__(self, parser, cwd, psFile):
        self.prsr = parser
        self.CWD  = cwd
        self.psFile = psFile
        self._notDiss = []
        self.internals = internals(self.prsr, self.CWD, self.psFile) 
        assert(not(hasattr(self.prsr, "internalName"))) 

    
    def getEquivalentBonds(self):
        if hasattr(self.prsr, "sampleSize"): 
            if self.prsr.AIMStype != "AIMS":
                tmpCWD = self.CWD + "/" + self.prsr.RNGdir + "1/" 
                tmpCWD += self.prsr.geomDir + "1" 
            else:
                tmpCWD = self.CWD + "/" + self.prsr.geomDir + "1" 
        else:
            tmpCWD = self.CWD  

        fileName = tmpCWD + "/Geometry.dat"
        dissPartners = self.prsr.dissPartners
        with open(fileName, "r") as geomLines:
            linenr = -2 
            atomNames = []
            for geomLine in geomLines:
                currAtmName = geomLine.strip().split()[0]
                linenr += 1
                if currAtmName in dissPartners:
                    atomNames.append((currAtmName,  str(linenr)))
                if linenr == 13:
                    break
            
            bondNames = [] 
            for atomName1 in atomNames:
                bondName = ""
                if atomName1[0] == dissPartners.split("-")[0]:
                    bondName += atomName1[0] + atomName1[1] + "-"
                    for atomName2 in atomNames:
                        if atomName2[0] == dissPartners.split("-")[1]:
                            bondNameN = bondName + atomName2[0] + atomName2[1]
                            bondNames.append(bondNameN)
            
            return bondNames


    def calculateAllInternals(self):
        """
            As described in the paper, the dissociation threshold
            is defined to be two times the standard deviation,
            assuming a normal distribution centred at the expected
            value. 
        """
        equivalentBonds = self.getEquivalentBonds()
        internal = []
        internalTime = []
        numSamples = 0
        if hasattr(self.prsr, "sampleSize"): 
            for geom in np.arange(1, self.prsr.sampleSize + 1):
                geomInternal  = []
                geomInternalT = []
                if geom in self.prsr.dupList:
                    continue
                if self.prsr.AIMStype != "AIMS":
                    for rng in np.arange(1, self.prsr.nrRNGs + 1):
                        rngInternal  = []
                        rngInternalT = []
                        tmpCWD = self.CWD + "/" + self.prsr.RNGdir + str(rng) 
                        tmpCWD += "/" + self.prsr.geomDir + str(geom)
                        for bond in equivalentBonds: 
                            equivInternal  = []
                            equivInternalT = []
                            self.internals.addInternal(equivInternal, tmpCWD, 
                                                       internalName = bond,
                                                       internalTime = equivInternalT)
                            rngInternal.append(equivInternal)
                            rngInternalT.append(equivInternalT)
                        geomInternal.append(rngInternal)
                        geomInternalT.append(rngInternalT)

                else:
                    tmpCWD = self.CWD + "/" + self.prsr.geomDir + str(geom)
                    for bond in equivalentBonds: 
                        equivInternal  = []
                        equivInternalT = []
                        self.internals.addInternal(equivInternal, tmpCWD, 
                                                   internalName = bond,
                                                   internalTime = equivInternalT)
                        geomInternal.append(equivInternal)
                        geomInternalT.append(equivInternalT)
                internal.append(geomInternal)
                internalTime.append(geomInternalT)
        else:
            tmpCWD = self.CWD 
            for bond in equivalentBonds: 
                equivInternal  = []
                equivInternalT = []
                self.internals.addInternal(equivInternal, tmpCWD,
                                           internalName = bond,
                                           internalTime = equivInternalT)
                self.addInternal(internal, tmpCWD)

        return internal, internalTime

    def addPopandBl(self, tmpTBFpops, geomTBFpop, geomBl, geomTime, internal,
                    internalTime):
        """
            Calculate means 
        """ 
        bondBlMean = []
        normalizer = 0
        maxDisp = 0
        maxBond = 0
        maxBondSpawn = []
        for i in np.arange(len(internal[0])):
            spawnMean = 0 
            for j in np.arange(len(internal)): 
                finDisp = internal[j][i][-1] - internal[j][i][0] 
                if finDisp > maxDisp:
                    maxDisp = finDisp
                    maxBond = j
            maxBondSpawn.append(maxBond)

        bondTBFpop = []
        bondTime = []
        bondBl = []
        for itmp, tmpTBFpop in enumerate(tmpTBFpops):
            tmpTime =  internalTime[maxBond][itmp]
            ntmpTime, ntmpTBFpop = self.psFile.zeroPadArray(tmpTime,
                                                            tmpTBFpop)
            interpTBFpop = np.interp(self.prsr.interpTime, ntmpTime,
                                     ntmpTBFpop)
            normalizer += interpTBFpop
        for iSpawn, spawn in enumerate(internal[maxBond]):
            blTBFpop = []
            blTime   = [] 
            blBl     = [] 
            bl0 = spawn[0] 
            for iBl, bl in enumerate(spawn):
                disp = bl - bl0 
                if (disp > (self.prsr.thresh)):
                    blBl.append(bl)
                    if tmpTBFpops[iSpawn].size == 1:
                        blTBFpop.append(tmpTBFpops[iSpawn])
                        blTime.append(internalTime[maxBond][iSpawn][0])
                    else:
                        blTBFpop.append(tmpTBFpops[iSpawn][iBl])
                        blTime.append(internalTime[maxBond][iSpawn][iBl])
            if len(blBl) > 0:
                blBl = np.array(blBl)
                blTime = np.array(blTime)
                blTBFpop = np.array(blTBFpop) 
                # discard those bonds that recross the dissociation threshold  
                if blTime[-1] != internalTime[maxBond][iSpawn][-1]:
                    continue
                nBlTime, blTBFpop = self.psFile.zeroPadArray(blTime,
                                                             blTBFpop)
                blTBFpop = (np.interp(self.prsr.interpTime, nBlTime,
                            blTBFpop) / normalizer)
                nBlTime, blBl = self.psFile.zeroPadArray(blTime, blBl)
                blBl = np.interp(self.prsr.interpTime, nBlTime, blBl)
                bondBl.append(blBl)
                bondTBFpop.append(blTBFpop)
                bondTime.append(nBlTime)
        totPop = 0
        for iPop, pop in enumerate(bondTBFpop):
            timeOld = 0
            for time in np.arange(pop.size):
                if bondBl[iPop][time] > 0.0:
                    timeNew =  self.prsr.interpTime[time]
                    diff = timeNew - timeOld
                    timeOld = timeNew
            totPop += np.array(pop)
        if isinstance(totPop, np.ndarray):
            for i in np.arange(totPop.size):
                if totPop[i] == 1.:
                    totPop[i:] = 1.
                    break 
            geomBl.append(bondBl)
            geomTBFpop.append(totPop)
            geomTime.append(bondTime)
        else:
            self._notDiss.append(-1)
        

    def getPopandBl(self, geomTBFpop, geomBl, geomTime, internal, internalTime): 
        TBFpops = self.psFile.getTBFpopulations(interp = False)
        ind = 0
        if hasattr(self.prsr, "sampleSize"):
            for geom in np.arange(1, self.prsr.sampleSize + 1):
                if geom in self.prsr.dupList:
                    continue
                if self.prsr.AIMStype != "AIMS":
                    rngTBFpop     = []
                    for rng in np.arange(self.prsr.nrRNGs):
                        tmpTBFpops    = TBFpops[ind][rng]
                        self.addPopandBl(tmpTBFpops, rngTBFpop, geomBl, 
                                         geomTime, internal[ind][rng], 
                                         internalTime[ind][rng])
                    geomTBFpop.append(rngTBFpop)
                    ind += 1
                else:
                    tmpTBFpops = TBFpops[ind]
                    self.addPopandBl(tmpTBFpops, geomTBFpop, geomBl, geomTime,
                                     internal[ind], internalTime[ind])
                    ind += 1


    def getGeomMolpop(self, geomMolpop, geomTBFpop, geomBl):
        if hasattr(self.prsr, "sampleSize"):
            for geom in np.arange(self.prsr.sampleSize-len(self.prsr.dupList)
                                  -len(self._notDiss)):
                if geom in self.prsr.dupList:
                    continue
                if self.prsr.AIMStype != "AIMS":
                    for rng in np.arange(len(geomTBFpop[geom])):
                        geomMolpop.append(geomTBFpop[geom][rng])
                else:
                    geomMolpop.append(geomTBFpop[geom])

    def globalExpec(self,internal,internalTime):
        equivalentBonds = self.getEquivalentBonds()
        numEquivBonds   = len(equivalentBonds)
        TBFpops = self.psFile.getTBFpopulations()
        TBFpopsN = self.psFile.getTBFpopulations(interp=False)
        samples = []
        for geom in np.arange(self.prsr.sampleSize-len(self.prsr.dupList)):
            for bond in np.arange(numEquivBonds):
                if self.prsr.AIMStype != "AIMS":
                    for rng in np.arange(self.prsr.nrRNGs):
                        tmpCWD = self.CWD + "/" + self.prsr.RNGdir + str(rng+1)  
                        tmpCWD += "/" + self.prsr.geomDir + str(geom+1)
                        tmpInternals = []
                        for spawn in np.arange(len(internal[geom][rng][bond])):
                            tmpInternal = np.interp(self.prsr.interpTime,  
                                    internalTime[geom][rng][bond][spawn],
                                    internal[geom][rng][bond][spawn])  
                            tmpInternals.append(tmpInternal)
                            self.internals.saveInternals(tmpCWD, internal[geom][rng][bond][spawn], spawn+1,
                                                         internalName = equivalentBonds[bond],
                                                         popTBF = TBFpopsN[geom][rng][spawn],
                                                         time = internalTime[geom][rng][bond][spawn])
                        sample = self.internals.calcIncoherentSum(
                                                    TBFpops[geom][rng],
                                                    tmpInternals)
                        samples.append(sample)
                else:
                    tmpInternals = []
                    tmpCWD = self.CWD + "/" + self.prsr.geomDir + str(geom+1)
                    for spawn in np.arange(len(internal[geom][bond])):
                        tmpInternal = np.interp(self.prsr.interpTime,  
                                     internalTime[geom][bond][spawn],
                                     internal[geom][bond][spawn])  
                        tmpInternals.append(tmpInternal)
                        #print len(TBFpopsN[geom])
                        self.internals.saveInternals(tmpCWD, internal[geom][bond][spawn], spawn+1,
                                                     internalName = equivalentBonds[bond],
                                                     popTBF = TBFpopsN[geom][spawn],
                                                     time = internalTime[geom][bond][spawn])
                    sample = self.internals.calcIncoherentSum(
                                                 TBFpops[geom],
                                                  tmpInternals)
                    samples.append(sample)

        
        mean = 0 
        for sample in samples:
            mean += sample
        mean = mean / len(samples)
        sd = 0 
        for sample in samples:
            sd += (sample - mean)**2 
        sd = np.sqrt(sd / (len(samples) * (len(samples) - 1)))
        saveFile  = self.prsr.internalType + self.prsr.dissPartners + "_"
        saveFile += self.prsr.AIMStype + ".dat"
        writeNPFile(3, saveFile, [self.prsr.interpTime, mean, sd], 
                    fmtStyle = "%8d %30.18e %30.18e")

    def getMolpop(self):
        internal, internalTime = self.calculateAllInternals()
        self.globalExpec(internal,internalTime)
        geomTBFpop = [] 
        geomBl     = []
        geomTime   = []
        self.getPopandBl(geomTBFpop, geomBl, geomTime, internal, internalTime)
        geomMolpop = []
        self.getGeomMolpop(geomMolpop, geomTBFpop, geomBl)
        mean = 0 
        for i in geomMolpop:
            mean += i
        totNum = len(geomMolpop) + len(self._notDiss)
        print(totNum)
        mean = mean / totNum
        sd = 0 
        for i in geomMolpop:
            sd += (i - mean)**2 
        sd = np.sqrt(sd / (totNum * (totNum - 1)))

        saveFile = "N_DISS_" + self.prsr.AIMStype + ".dat"
        writeNPFile(3, saveFile, [self.prsr.interpTime, mean, sd], 
                    fmtStyle = "%8d %30.18e %30.18e")
        saveFile = "N_UNDISS_" + self.prsr.AIMStype + ".dat"
        writeNPFile(3, saveFile, [self.prsr.interpTime, 1 - mean, sd], 
                    fmtStyle = "%8d %30.18e %30.18e")
