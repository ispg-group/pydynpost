#!/usr/bin/env python
import numpy as np
import os
import math
from filesys import *
from commonmethods.misc import *
from commonmethods.parse import *
from commonmethods.internals import *
from aimsscripts.aimssave import *

class internals(object): 
    def __init__(self, parser, cwd, psFile):
        self.prsr = parser
        self.CWD  = cwd
        self.psFile = psFile

    def addInternal(self, positions, internal, internalTime,
                    internalName = None):
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

        for molStruct in positions[0]:
            atmcoords = [] 
            for atm in molStruct[1]:
                for i in np.arange(len(atmsInvolved)):
                    if atm[0] == atmsInvolved[i]:
                        atmcoords.append(np.array(atm[1:]))
            internal.append(calcIntrnl(atmcoords, self.prsr.internalType))
            internalTime.append(molStruct[0])

    
    def getInternals(self):
        positions = self.psFile.readPositions("movie.xyz", "geom")


class molpop(object): 
    def __init__(self, parser, cwd, psFile):
        self.prsr = parser
        self.CWD  = cwd
        self.psFile = psFile
        self.internals = internals(self.prsr, self.CWD, self.psFile)

    def getEquivalentBonds(self):
        if not(self.prsr.nrRNGs == 0):
            tmpCWD = self.CWD + "/" + self.prsr.RNGdir + "1/" 
            tmpCWD += self.prsr.geomDir + "1" 
        else:
            tmpCWD = self.CWD + "/" + self.prsr.geomDir + "1" 

        fileName = tmpCWD + "/geom"
        dissPartners = self.prsr.dissPartners
        with open(fileName, "r") as geomLines:
            linenr = 0 
            atomNames = []
            geomLines.readline()    
            geomLines.readline()    
            for geomLine in geomLines:
                currAtmName = geomLine.strip().split()[0]
                linenr += 1
                if currAtmName in dissPartners:
                    atomNames.append((currAtmName,  str(linenr)))
            
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
        positions = self.psFile.readPositions("movie.xyz", "geom")
        #print len(positions[0])
        #print hasattr(self, "addInternal")
        internal = []
        internalTime = []
        numSamples = 0
        if hasattr(self.prsr, "sampleSize"): 
            for geom in np.arange(self.prsr.sampleSize
                                  -len(self.prsr.dupList)):
                geomInternal  = []
                geomInternalT = []
                if not(self.prsr.nrRNGs == 0):
                    for rng in np.arange(self.prsr.nrRNGs):
                        rngInternal  = []
                        rngInternalT = []
                        for bond in equivalentBonds: 
                            equivInternal  = []
                            equivInternalT = []
                            self.internals.addInternal(positions[geom][rng], equivInternal,  
                                                       equivInternalT, internalName = bond)
                            saveFile = self.CWD + "/" + self.prsr.RNGdir + str(rng + 1) 
                            saveFile += "/" + self.prsr.geomDir + str(geom + 1) + "/"
                            saveFile += self.prsr.internalType + bond + ".dat" 
                            writeNPFile(2, saveFile, [np.array(equivInternalT)*0.02418884254, 
                                                      np.array(equivInternal)], 
                                        fmtStyle = "%8.2f %30.18e")
                            rngInternal.append(equivInternal)
                            rngInternalT.append(equivInternalT)
                        geomInternal.append(rngInternal)
                        geomInternalT.append(rngInternalT)

                else:
                    tmpCWD = self.CWD + "/" + self.prsr.geomDir + str(geom)
                    for bond in equivalentBonds: 
                        equivInternal  = []
                        equivInternalT = []
                        self.internals.addInternal(positions[geom], equivInternal,  
                                                   equivInternalT, internalName = bond)
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

    def calculateMolpop(self, internal, internalTime):
        numDissoc = np.zeros(self.prsr.interpTime.size) 
        for geom in np.arange(self.prsr.sampleSize
                              -len(self.prsr.dupList)):
            if not(self.prsr.nrRNGs == 0):
                for rng in np.arange(self.prsr.nrRNGs):
                    maxDisp = 0 
                    maxBond = 0
                    for bond in np.arange(len(internal[geom][rng])):
                        tmpInternal = internal[geom][rng][bond] 
                        tmpInternalTime = internalTime[geom][rng][bond] 
                        tmpDissoc = np.zeros(len(tmpInternal)) 
                        blTime = []
                        blBl = []
                        bl0 = tmpInternal[0]
                        for iBl, bl in enumerate(tmpInternal): 
                            disp = bl - bl0 
                            if disp > maxDisp:
                                maxBond = bond 
                                maxDisp = disp

                    maxInternal = internal[geom][rng][maxBond] 
                    maxInternalTime = internalTime[geom][rng][maxBond] 
                    tmpDissoc = np.zeros(len(tmpInternal)) 
                    blTime = []
                    blBl = []
                    bl0 = maxInternal[0]
                    for iBl, bl in enumerate(maxInternal): 
                        disp = bl - bl0 
                        if (disp >  self.prsr.thresh):
                            tmpDissoc[iBl:] += 1
                            blTime.extend(tmpInternalTime[iBl:])
                            break
                    
                    # discard those bonds that recross the dissociation threshold  
                    if len(blTime) != 0:
                        if blTime[-1] != tmpInternalTime[-1]:
                            continue
                        else:
                            told = 0.0 
                            iTold = 0
                        tmpDissoc = np.interp(self.prsr.interpTime, 
                                              tmpInternalTime, tmpDissoc)
                        numDissoc += tmpDissoc

        numSamples = self.prsr.nrRNGs * (self.prsr.sampleSize - len(self.prsr.dupList))
        mean = numDissoc / numSamples
        dissFile = "N_DISS.dat" 
        np.savetxt(dissFile, np.array([self.prsr.interpTime, mean]).T,
                   fmt="%8d %30.18e")
        undissFile = "N_UNDISS.dat" 
        np.savetxt(undissFile, np.array([self.prsr.interpTime, 1-mean]).T,
                   fmt="%8d %30.18e")
                                

    def getMolpop(self):
        internal, internalTime = self.calculateDissociationThresh()
        self.calculateMolpop(internal, internalTime)
