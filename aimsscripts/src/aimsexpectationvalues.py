#!/usr/local/Cluster-Apps/python/2.7.9/bin/python
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

class observables(object):
    """
        This abstract base class calculates different observables 
        as incoherent sums (as of 10/2021) over TBF populations. 
        A derived class just has to implement the getRawObservable
        method, and call the template method getExpectationValue
        from its own method called getObservable (where Observable
        can be Internal, DipoleMoment, etc.)
   """
    __metaclass__ = ABCMeta
    def __init__(self, parser, cwd, psFile):
        self.prsr = parser
        self.CWD  = cwd
        self.psFile = psFile

    def calcIncoherentSum(self, TBFpops, observable): 
        # It simply do what it do!
        incoherentSum = np.zeros(TBFpops[0].size)
        numerator     = np.zeros(TBFpops[0].size)
        denominator   = np.zeros(TBFpops[0].size)
        for i in np.arange(len(TBFpops)):
            tmpTerm = TBFpops[i] * observable[i]
            numerator += tmpTerm
            denominator += TBFpops[i]

        incoherentSum = numerator / denominator
        return incoherentSum

    @abc.abstractmethod
    def getRAWObservable(self): 
        pass

    def getExpectationValue(self, saveStringParams, internalName = None):
        expecObservable = np.zeros(self.prsr.interpTime.size) 
        TBFpops   = self.psFile.getTBFpopulations()
        rawObservable = self.getRAWObservable()
        geomExpects = []
        ind = 0 
        if hasattr(self.prsr, "sampleSize"):
            for geom in np.arange(1, self.prsr.sampleSize + 1):
                if geom in self.prsr.dupList:
                    #print geom
                    continue
                if self.prsr.AIMStype != "AIMS":
                    geomExpecObservable = np.zeros(self.prsr.interpTime.size)
                    tmpExpects = []  
                    for rng in np.arange(self.prsr.nrRNGs):
                        tmpTBFpops    = TBFpops[ind][rng]
                        tmpObservable = rawObservable[ind][rng]
                        tmpExpect     = self.calcIncoherentSum(tmpTBFpops, 
                                                               tmpObservable)
                        tmpExpects.append(tmpExpect)
                        geomExpects.append(tmpExpect)
                        expecObservable += tmpExpect
                        geomExpecObservable += tmpExpect
                        #tmpCWD  = self.CWD + "/rng" + str(rng + 1) + "/geom_" 
                        tmpCWD  = self.CWD + "/" + self.prsr.RNGdir + str(rng + 1)  
                        #tmpCWD += str(geom + 1)
                        tmpCWD += "/" + self.prsr.geomDir + str(geom + 1)
                        saveString = ""
                        for saveStringParam in saveStringParams: 
                            saveString += saveStringParam 
                        if not(internalName == None):
                            saveString += internalName[ind]
                        saveString += ".dat" 
                        saveFile   = self.psFile.inputFileName(tmpCWD, saveString) 
                        writeNPFile(2, saveFile, [self.prsr.interpTime, tmpExpect],
                                    fmtStyle = "%8d %30.18e")

                    geomExpecObservable = geomExpecObservable / self.prsr.nrRNGs
                    geomStdErrObservable = np.zeros(geomExpecObservable.size)
                    for i in np.arange(len(tmpExpects)):
                        geomStdErrObservable += (tmpExpects[i] 
                                                 - geomExpecObservable)**2
                    geomStdErrObservable = np.sqrt(geomStdErrObservable
                                                 / (self.prsr.nrRNGs 
                                                    * (self.prsr.nrRNGs + 1)))
                    tmpCWD = self.CWD 
                    saveString = "" 
                    for saveStringParam in saveStringParams:
                        saveString += saveStringParam
                    saveString += "_IC" + str(geom + 1) + "_" + str(rng + 1)
                    saveString += ".dat"
                    saveFile   = self.psFile.inputFileName(tmpCWD, saveString) 
                    writeNPFile(3, saveFile, [self.prsr.interpTime, 
                                              geomExpecObservable,
                                              geomStdErrObservable], 
                                fmtStyle = "%8d %30.18e %30.18e")
                    ind += 1
                else:
                    tmpTBFpops    = TBFpops[ind]
                    tmpObservable = rawObservable[ind]
                    ind += 1
                    tmpExpect     = self.calcIncoherentSum(tmpTBFpops, 
                                                           tmpObservable)
                    expecObservable += tmpExpect
                    geomExpects.append(tmpExpect)
                    #tmpCWD = self.CWD + "/geom_" + str(geom + 1)
                    tmpCWD = self.CWD + "/" + self.prsr.geomDir + str(geom)
                    saveString  = "" 
                    for saveStringParam in saveStringParams:
                        saveString += saveStringParam
                    if not(internalName == None):
                        saveString += internalName[ind]
                    saveString += ".dat"
                    #print tmpCWD, saveString
                    saveFile    = self.psFile.inputFileName(tmpCWD, saveString)
                    writeNPFile(2, saveFile, [self.prsr.interpTime, tmpExpect],
                                fmtStyle = "%8d %30.18e")

            if self.prsr.AIMStype == "AIMS":
                Nsamples = self.prsr.sampleSize - len(self.prsr.dupList)
                if hasattr(self.prsr, "internalName"):
                    saveFile  = self.CWD + "/" + self.prsr.internalType 
                    saveFile += self.prsr.internalName 
                else:
                    saveFile = self.CWD + "/" + self.prsr.internalType 
                saveFile += "_AIMS.dat"
            else:
                Nsamples  = (self.prsr.sampleSize - len(self.prsr.dupList)) 
                Nsamples *= self.prsr.nrRNGs
                if hasattr(self.prsr, "internalName"):
                    saveFile  = self.CWD + "/" + self.prsr.internalType
                    saveFile += self.prsr.internalName
                else:
                    saveFile = self.CWD + "/" + self.prsr.internalType
                saveFile += "_" + self.prsr.AIMStype + "_" 
                saveFile += str(self.prsr.nrRNGs) + ".dat"
            expecObservable = expecObservable / Nsamples
            stdErrObservable = np.zeros(expecObservable.size)
            for i in np.arange(len(geomExpects)):
                stdErrObservable += (geomExpects[i] - expecObservable)**2
            stdErrObservable = np.sqrt(stdErrObservable / (Nsamples * 
                               (Nsamples + 1)))
            writeNPFile(3, saveFile, [self.prsr.interpTime, expecObservable,
                                      stdErrObservable], 
                        fmtStyle = "%8d %30.18e %30.18e")
        else:
            expecObservable  = self.calcIncoherentSum(TBFpops, rawObservable)
            tmpCWD = self.CWD 
            saveString  = "" 
            for saveStringParam in saveStringParams:
                saveString += saveStringParam
            saveString += ".dat"
            saveFile    = self.psFile.inputFileName(tmpCWD, saveString) 
            writeNPFile(2, saveFile, [self.prsr.interpTime, expecObservable],
                        fmtStyle = "%8d %30.18e")

class internals(observables):
    def __init__(self, parser, cwd, psFile, internalName = None):
        observables.__init__(self, parser, cwd, psFile)
        self.observableType = "internals" 
        self.internalName = internalName

    def getInternals(self):
        internals = 0
        if self.internalName == None:
            internals = self.getExpectationValue(saveStringParams = 
                                                 [self.prsr.internalType,
                                                  self.prsr.internalName])
        else:
            internals = self.getExpectationValue(saveStringParams = 
                                                 [self.prsr.internalType],
                                                 internalName = 
                                                 self.internalName) 

    def saveInternals(self, tmpCWD, curInternals, ID, internalName = None,
                      popTBF = None, time = None):
        if internalName == None:
            saveString  = self.prsr.internalType + self.prsr.internalName
        else: 
            saveString  = self.prsr.internalType + internalName
        saveString += "." + str(ID) 
        saveFile    = self.psFile.inputFileName(tmpCWD, saveString)
         
        if ((popTBF == None) and (time == None)): 
            saveInternals   = curInternals[curInternals > 0]
            saveTime        = self.prsr.interpTime[curInternals > 0]
            writeNPFile(2, saveFile, [saveTime, saveInternals], 
                        fmtStyle = "%8d %30.18e") 
        else:
            #print popTBF
            if not(len(time) == 1):
                saveTime        = np.linspace(time[0], time[-1], 1000)
                saveInternals   = np.interp(saveTime, time, curInternals)
                savePop         = np.interp(saveTime, time, popTBF)
                writeNPFile(3, saveFile, [saveTime, saveInternals, savePop],
                            fmtStyle = "%8d %30.18e %30.18e") 

    def calcInternals(self, coords, internalName = None):
        if internalName == None:
            atmsInvolved = self.prsr.internalName.split("-") 
        else: 
            atmsInvolved = internalName.split("-") 
        if self.prsr.internalType == "bl":
            assert (len(atmsInvolved) == 2)
        elif self.prsr.internalType == "ba":
            assert (len(atmsInvolved) == 3)
        elif self.prsr.internalType == "td":
            assert (len(atmsInvolved) == 4)
        intrnl     = []
        intrnlTime = []
        for molStruct in coords: 
            atmcoords = [] 
            for atm in molStruct[1]:
                for i in np.arange(len(atmsInvolved)):
                    if atm[0] == atmsInvolved[i]:
                        atmcoords.append(np.array(atm[1:]))
            intrnl.append(calcIntrnl(atmcoords, self.prsr.internalType))
            intrnlTime.append(molStruct[0])
        intrnl     = np.array(intrnl)
        intrnlTime = np.array(intrnlTime)
        intrnlTime, intrnl = self.psFile.zeroPadArray(intrnlTime, intrnl)

        return intrnl, intrnlTime

    def addInternal(self, internal, tmpCWD, ind = None, internalName = None,
                    internalTime = None):
        numParticles = self.psFile.findNumAtoms(tmpCWD)
        spawnTimes, childIDs, parentIDs, numSpawns = self.psFile.findNrSpawns(tmpCWD)
        posErr, FGcoord = self.psFile.readPositions(1, tmpCWD,
                                                    numParticles)
        if ((self.internalName == None) and 
           (internalName == None)):
            FGinternal, FGtime = self.calcInternals(FGcoord)
        elif ((self.internalName == None) and 
              (internalName != None)):
            FGinternal, FGtime = self.calcInternals(FGcoord, 
                                                    internalName =
                                                    internalName)
        else:
            FGinternal, FGtime = self.calcInternals(FGcoord, 
                                 internalName = self.internalName[ind])
        interpFGintrnl = np.interp(self.prsr.interpTime, FGtime, 
                                   FGinternal)
        if (internalName == None):
            internal.append(interpFGintrnl)
        else:
            internalTime.append(FGtime[FGinternal > 0])
            internal.append(FGinternal[FGinternal > 0])
        if ((self.internalName == None) and
            (internalName == None)):
            self.saveInternals(tmpCWD, interpFGintrnl, 1)
       # elif ((self.internalName == None) and
       #     (internalName != None)):
       #     self.saveInternals(tmpCWD, interpFGintrnl, 1,
       #                        internalName = internalName)
        elif ((self.internalName != None) and
              (internalName == None)):
            self.saveInternals(tmpCWD, interpFGintrnl, 1, 
                               internalName = self.internalName[ind])
         
        if type(childIDs) != np.ndarray:
            childIDs = np.array([childIDs])
        for childID in childIDs:
            posErr, CHcoord = self.psFile.readPositions(childID, 
                                                        tmpCWD,
                                                        numParticles)
            if not(posErr):
                if ((self.internalName == None) and 
                   (internalName == None)):
                    CHinternal, CHtime = self.calcInternals(CHcoord)
                elif ((self.internalName == None) and 
                      (internalName != None)):
                    CHinternal, CHtime = self.calcInternals(CHcoord, 
                                                            internalName =
                                                            internalName)
                else:
                    CHinternal, CHtime = self.calcInternals(CHcoord,
                                         internalName = self.internalName[ind])
                interpCHintrnl = np.interp(self.prsr.interpTime,  
                                           CHtime, CHinternal)
                if ((self.internalName == None) and
                    (internalName == None)):
                        self.saveInternals(tmpCWD, interpCHintrnl, childID)
       #         elif ((self.internalName == None) and
       #             (internalName != None)):
       #                 self.saveInternals(tmpCWD, interpCHintrnl, childID,
       #                                    internalName = internalName)
                elif ((self.internalName != None) and
                      (internalName == None)):
                    self.saveInternals(tmpCWD, interpFGintrnl, childID, 
                                       internalName = self.internalName[ind])

                if (internalName == None):
                    internal.append(interpCHintrnl)
                else:
                    internalTime.append(CHtime[CHinternal > 0])
                    internal.append(CHinternal[CHinternal > 0])

    def getRAWObservable(self): 
        internal = []
        if hasattr(self.prsr, "sampleSize"): 
            if not(self.internalName == None):
                ind = 0
            for geom in np.arange(1, self.prsr.sampleSize + 1):
                geomInternal = []
                if geom in self.prsr.dupList:
                    #print geom
                    continue
                if self.prsr.AIMStype != "AIMS":
                    for rng in np.arange(1, self.prsr.nrRNGs + 1):
                        rngInternal = []
                        tmpCWD = self.CWD + "/" + self.prsr.RNGdir + str(rng)  
                        tmpCWD += "/" + self.prsr.geomDir + str(geom)
                        if (self.internalName == None):
                            self.addInternal(rngInternal, tmpCWD)
                        else:
                            self.addInternal(rngInternal, tmpCWD, ind = ind)
                        geomInternal.append(rngInternal)

                    if not(self.internalName == None):
                        ind += 1
                else:
                    #tmpCWD = self.CWD + "/geom_" + str(geom)
                    tmpCWD = self.CWD + "/" + self.prsr.geomDir + str(geom)
                    if (self.internalName == None):
                        self.addInternal(geomInternal, tmpCWD) 
                    else:
                        self.addInternal(geomInternal, tmpCWD, ind = ind) 
                if len(geomInternal) != 0:
                    internal.append(geomInternal)
        else:
            tmpCWD = self.CWD 
            self.addInternal(internal, tmpCWD)

        return internal

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
        self.internals = internals(self.prsr, self.CWD, self.psFile) 
        assert(not(hasattr(self.prsr, "internalName"))) 

    #def getDissociation(self):
    #    self.internals.getInternals()
    
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
        #numParticles = self.psFile.findNumAtoms(tmpCWD)
        #test = redundantInternals(numParticles, fname = fileName,
        #                          fileType = "dat")
        #test.getBondPartner()
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


    def calculateDissociationThresh(self):
        """
            As described in the paper, the dissociation threshold
            is defined to be two times the standard deviation,
            assuming a normal distribution centred at the expected
            value. 
        """
        equivalentBonds = self.getEquivalentBonds()
        #print hasattr(self, "addInternal")
        internal = []
        internalTime = []
        numSamples = 0
        if hasattr(self.prsr, "sampleSize"): 
            for geom in np.arange(1, self.prsr.sampleSize + 1):
                geomInternal  = []
                geomInternalT = []
                if geom in self.prsr.dupList:
                    #print geom
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
        #print(len(internal[0]), len(internal))
        maxDisp = 0
        maxBond = 0
        maxBondSpawn = []
        for i in np.arange(len(internal[0])):
            spawnMean = 0 
            for j in np.arange(len(internal)): 
                currIntrnl = np.array(internal[j][i]) - internal[j][i][0] 
                for k in currIntrnl:
                    if k > maxDisp:
                        maxDisp = k
                        maxBond = j
                spawnMean += currIntrnl 
            maxBondSpawn.append(maxBond)
            #print maxBondSpawn
            spawnMean = spawnMean / len(internal)
            bondBlMean.append(spawnMean)

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
                #if (disp > (bondBlMean[iSpawn][iBl] + self.prsr.thresh)):
                #print disp, self.prsr.thresh, internalTime[maxBond][iSpawn][iBl], maxBond + 1, iSpawn + 1
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
                    if  timeOld > 0.0 and diff > self.prsr.step: 
                        print diff, timeOld, timeNew
                    timeOld = timeNew
            totPop += np.array(pop)
        for i in np.arange(totPop.size):
            if totPop[i] == 1.:
                totPop[i:] = 1.
                break 
        geomBl.append(bondBl)
        geomTBFpop.append(totPop)
        geomTime.append(bondTime)
        #print geomTime
        

    def getPopandBl(self, geomTBFpop, geomBl, geomTime, internal, internalTime): 
        TBFpops = self.psFile.getTBFpopulations(interp = False)
        ind = 0
        if hasattr(self.prsr, "sampleSize"):
            for geom in np.arange(1, self.prsr.sampleSize + 1):
                if geom in self.prsr.dupList:
                    #print geom
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
                    #print(geom)
                    ind += 1
            #print len(geomTBFpop)

    #def addGeomMolpop(self, geomMolpop, geomTBFpop):
    #    if len(geomTBFpop) == 0:
    #        geomMolpop.append(np.zeros(self.prsr.interpTime.size))   
    #    else:
    #        bondMolpops = []
    #        numSamples = 0
    #        for iBond, bond in enumerate(geomTBFpop):
    #            #print bond
    #            tmpBondMolpop = 0 
    #            for iSpawn, spawn in enumerate(bond):
    #                tmpBondMolpop += spawn  
    #            bondMolpops.append(tmpBondMolpop) 
    #        if len(bondMolpops) > 1:
    #            #print "more than one"
    #            #print bondMolpops
    #            bondMolpop = np.zeros(self.prsr.interpTime.size) 
    #            for itime in np.arange(self.prsr.interpTime.size):
    #                average = [] 
    #                for bond in bondMolpops: 
    #                    if bond[itime] > 0.0:
    #                        average.append(True)
    #                    else:
    #                        average.append(False)
    #                denominator = sum(average)
    #                for bond in bondMolpops:
    #                    if denominator == 0:
    #                        bondMolpop[itime] += bond[itime]
    #                    #bondMolpop[itime] += bond[itime]
    #                    else:
    #                        bondMolpop[itime] += bond[itime]/denominator
    #                        
    #        else:
    #            #print "just one"
    #            bondMolpop = tmpBondMolpop 
    #        geomMolpop.append(bondMolpop)

    def getGeomMolpop(self, geomMolpop, geomTBFpop, geomBl):
        #print(len(geomTBFpop))
        if hasattr(self.prsr, "sampleSize"):
            for geom in np.arange(self.prsr.sampleSize):
                if geom in self.prsr.dupList:
                    #print geom
                    continue
                if self.prsr.AIMStype != "AIMS":
                    for rng in np.arange(self.prsr.nrRNGs):
                        #print len(geomTBFpop[geom][rng])
                        geomMolpop.append(geomTBFpop[geom][rng])
                else:
                    geomMolpop.append(geomTBFpop[geom])

    def globalExpec(self,internal,internalTime):
        equivalentBonds = self.getEquivalentBonds()
        numEquivBonds   = len(equivalentBonds)
        #print numEquivBonds
        TBFpops = self.psFile.getTBFpopulations()
        TBFpopsN = self.psFile.getTBFpopulations(interp=False)
        samples = []
        for geom in np.arange(self.prsr.sampleSize):
            for bond in np.arange(numEquivBonds):
                if self.prsr.AIMStype != "AIMS":
                    for rng in np.arange(self.prsr.nrRNGs):
                        tmpCWD = self.CWD + "/" + self.prsr.RNGdir + str(rng+1)  
                        tmpCWD += "/" + self.prsr.geomDir + str(geom+1)
                        tmpInternals = []
                        for spawn in np.arange(len(internal[geom][rng][bond])):
                            #if (spawn + 1) == len(internal[geom][rng][bond]):
                            #    print geom+1, rng+1, spawn+1
                            #    print TBFpops[geom][rng] 
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
        #print len(geomMolpop)
        sd = 0 
        for sample in samples:
            sd += (sample - mean)**2 
        sd = np.sqrt(sd / (len(samples) * (len(samples) - 1)))
        saveFile  = self.prsr.internalType + self.prsr.dissPartners + "_"
        saveFile += self.prsr.AIMStype + ".dat"
        writeNPFile(3, saveFile, [self.prsr.interpTime, mean, sd], 
                    fmtStyle = "%8d %30.18e %30.18e")

    def getMolpop(self):
        internal, internalTime = self.calculateDissociationThresh()
        self.globalExpec(internal,internalTime)
        geomTBFpop = [] 
        geomBl     = []
        geomTime   = []
        self.getPopandBl(geomTBFpop, geomBl, geomTime, internal, internalTime)
        geomMolpop = []
        self.getGeomMolpop(geomMolpop, geomTBFpop, geomBl)
        mean = 0 
        for i in geomMolpop:
            #print i
            mean += i
        mean = mean / len(geomMolpop)
        sd = 0 
        for i in geomMolpop:
            sd += (i - mean)**2 
        sd = np.sqrt(sd / (len(geomMolpop) * (len(geomMolpop) - 1)))

        saveFile = "N_DISS_" + self.prsr.AIMStype + ".dat"
        writeNPFile(3, saveFile, [self.prsr.interpTime, mean, sd], 
                    fmtStyle = "%8d %30.18e %30.18e")
        saveFile = "N_UNDISS_" + self.prsr.AIMStype + ".dat"
        writeNPFile(3, saveFile, [self.prsr.interpTime, 1 - mean, sd], 
                    fmtStyle = "%8d %30.18e %30.18e")
