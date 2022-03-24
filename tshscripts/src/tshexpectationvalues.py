#!/usr/local/Cluster-Apps/python/2.7.9/bin/python
import numpy as np
import os
import math
from filesys import *
from misc import *
from parse import *
from internals import *

class internals(object): 
    def __init__(self, parser, cwd, psFile):
        self.prsr = parser
        self.CWD  = cwd
        self.psFile = psFile

    def addInternals(self, positions, internalName = None):
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

        internal     = []
        internalTime = []
        for geom in np.arange(self.prsr.sampleSize):
            geomIntrnl = []
            geomIntrnlTime = []
            if not(self.prsr.nrRNGs == 0):
                for rng in np.arange(self.prsr.sampleSize): 
                    rngIntrnl     = []
                    rngIntrnlTime = []
                    for molStruct in positions[geom][rng][0]:
                        atmcoords = [] 
                        for atm in molStruct[1]:
                            for i in np.arange(len(atmsInvolved)):
                                if atm[0] == atmsInvolved[i]:
                                    atmcoords.append(np.array(atm[1:]))
                        rngIntrnl.append(calcIntrnl(atmcoords, self.prsr.internalType))
                        rngIntrnlTime.append(molStruct[0])
                geomIntrnl.append(rngIntrnl)
                geomIntrnlTime.append(rngIntrnlTime)
            else:
                for molstruct in positions[geom][0]:
                    atmcoords = [] 
                    for atm in molStruct[1]:
                        for i in np.arange(len(atmsInvolved)):
                            if atm[0] == atmsInvolved[i]:
                                atmcoords.append(np.array(atm[1:]))
                    geomIntrnl.append(calcIntrnl(atmcoords, self.prsr.internalType))
                    geomIntrnlTime.append(molStruct[0])
            internal.append(geomIntrnl)
            internalTime.append(geomIntrnlTime)

        return internal, internalTime
    
    def getInternals(self):
        positions = self.psFile.readPositions()
        internals = self.addInternals(positions)


class molpop(object): 
    def __init__(self, parser, cwd, psFile):
        self.prsr = parser
        self.CWD  = cwd
        self.psFile = psFile
        self.intrnls = internals(self.prsr, self.CWD, self.psFile)

    def getMolpop(self):
