#!/usr/bin/env python
import numpy as np
import os
from commonmethods.filesys import *
from commonmethods.misc import *
from commonmethods.parse import *
from .aimsinp import *
from commonmethods.writefiles import *


class analyzeCouplings(object):
    """ 
        Class handling the analysis of effective
        nonadiabatic coupling parameters. Commonly
        used to judge the adequacy of the coupling
        threshold. 
    """
    def __init__(self, parser, cwd, psFile):
        self.prsr = parser
        self.CWD  = cwd
        self.psFile = psFile

    def saveCouplings(self, ID, couplings, tmpCWD):
        """ 
            Save the normal couplings to their 
            respective directory (especially important
            when coupV coupling type is used)
        """
        content = [couplings[0]]
        for i in couplings[1].T:
            content.append(i) 
        fileName = tmpCWD + "/coup." + str(ID)
        writeNPFile(self.prsr.nrStates + 1, fileName, content)

    def saveModCouplings(self, geom, ID, couplings, tmpCWD):
        """ 
            Save the modified couplings to their 
            respective directory
        """
        printContent = []
        fmt="%8d %30.18e"
        for i in np.arange(len(couplings[0])):
            printContent = [couplings[0][i], couplings[1][i]]
            fileName = tmpCWD + "/coupMod." + str(geom) + "."  
            fileName += str(ID) + "." + str(i + 1) 
            writeNPFile(2, fileName, printContent, fmtStyle = fmt)
            

    def saveAvCouplingProfile(self, avgCoupling):
        """ 
            Save the average coupling profile to
            the current working directory. 
        """
        fmt="%8d %30.18e"
        fileName = self.CWD + "/coupAvg.dat" 
        writeNPFile(2, fileName, avgCoupling, fmtStyle = fmt)

    def modifyCouplings(self, rCouplings, CSThresh):
        coupTime = rCouplings[0]  
        coupRaw  = rCouplings[1]
        mErr     = False
        modCTimes = []
        modCoups = []
        for cr in coupRaw.T: 
            if not(np.all(cr == 0)): 
                mask = (cr > CSThresh) 
            else:
                continue
            
            if np.any(mask):
                start  = False 
                starti = 0
                end    = True
                endi   = 0
                for i, mvl in enumerate(mask):
                    if (mvl == True and not(start)):
                        end   = False
                        start = True
                        starti = i
                        continue

                    if (mvl == False and not(end)):
                        endi  = i
                        end   = True
                        start = False
                        modCoup  = cr[starti:endi]
                        tmpTime  = coupTime[starti:endi]
                        imaxCoup = np.argmax(modCoup)
                        modCTime = tmpTime - tmpTime[imaxCoup] 
                        modCTimes.append(modCTime)
                        modCoups.append(modCoup)
            else:
                continue


        if ((len(modCTimes) == 0) or (len(modCoups) == 0)):
            mErr = True
        return mErr, [modCTimes, modCoups]  

    def getModCouplings(self): 
        """ 
            Function collates "modified" couplings, 
            i.e., those that are above the coupling 
            threshold and shifts their entry time 
            to zero.
        """
        modCouplings   = []
        for geom in np.arange(1, self.prsr.sampleSize + 1):
            if geom in self.prsr.dupList:
                continue
            geomModCouplings = []
            pCouplings = []
            tmpCWD = self.CWD + "/" + self.prsr.geomDir + str(geom)
            if not(hasattr(self.prsr, "CSThresh")):
                mErr, CSThresh = self.psFile.readCSThresh(tmpCWD) 
                assert(not(mErr))
            else:
                CSThresh = self.prsr.CSThresh
            spawnTimes, childIDs, parentIDs, numSpawns = self.psFile.findNrSpawns(tmpCWD)
            mErr, pCouplings = self.psFile.readCouplings(1, tmpCWD)
            if mErr:
                continue   
            mErr, tmpModCouplings = self.modifyCouplings(pCouplings, CSThresh)
            if mErr:
                continue   
            self.saveCouplings(1, pCouplings, tmpCWD)
            self.saveModCouplings(geom, 1, tmpModCouplings, tmpCWD)
            geomModCouplings.append(tmpModCouplings)
            if (childIDs.size != 1):
                for childID in childIDs:
                    mErr, cCouplings = self.psFile.readCouplings(childID, tmpCWD)
                    if mErr:
                        continue   
                    mErr, tmpModCouplings = self.modifyCouplings(cCouplings, CSThresh)
                    if mErr:
                        continue   
                    geomModCouplings.append(tmpModCouplings)
                    self.saveCouplings(childID, pCouplings, tmpCWD)
                    self.saveModCouplings(geom, childID, tmpModCouplings, tmpCWD)
                    cCouplings = []
            else:
                mErr, cCouplings = self.psFile.readCouplings(childIDs, tmpCWD)
                if mErr:
                    continue   
                mErr, tmpModCouplings = self.modifyCouplings(cCouplings, CSThresh)
                if mErr:
                    continue   
                geomModCouplings.append(tmpModCouplings)
                self.saveCouplings(childIDs, pCouplings, tmpCWD)
                self.saveModCouplings(geom, childIDs, tmpModCouplings, tmpCWD)
                cCouplings = []
            modCouplings.append(geomModCouplings)
        return modCouplings
            

    def getCoupling(self):
        """
            Main function gathering the modified couplings
            and returning the average coupling profile  
        """
        couplings = self.getModCouplings()
        # average the modified couplings
        infTime = 0
        supTime = 0
        for modCoupling in couplings:
            for geomCoupling in modCoupling:
                for coupling in geomCoupling[0]:
                    if coupling[0] <= infTime:
                        infTime = coupling[0]
                        #print infTime
                    if coupling[-1] >= supTime:
                        supTime = coupling[-1]
                        #print supTime
        nSteps = int((supTime - infTime)/self.prsr.step)
        tInterp = np.linspace(infTime, supTime, num = nSteps) 
        avgCoupling = np.zeros(tInterp.size) 
        numSamp     = 0
        for modCoupling in couplings:
            for geomCoupling in modCoupling:
                for i, coupling in enumerate(geomCoupling[1]):
                    cInterp = np.interp(tInterp, geomCoupling[0][i], coupling)  
                    avgCoupling += cInterp 
                    numSamp     += 1
        avgCoupling = avgCoupling / numSamp 
        # then save them in some file 
        self.saveAvCouplingProfile([tInterp, avgCoupling])

