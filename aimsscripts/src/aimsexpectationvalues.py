#!/usr/bin/env python
import numpy as np
import abc
from abc import ABCMeta
import os
import copy
from filesys import *
from misc import *
from parse import *
from aimsinp import *
from writefiles import *
from internals import *
from scipy import integrate
b2A = 0.529177249
A2b = 1./b2A

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
        """
             It simply do what it do!
        """
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
        """ 
            This function has to be implemented by all derived 
            classes. The result is needed to calculate the 
            incoherent sum step of the general getExpectationValue 
            method.
        """
        pass

    @abc.abstractmethod
    def calcRAWObservable(self, molStruct):
        """ 
            This function has to be implemented by all derived 
            classes. The result is needed to calculate the coherent
            sum via monte carlo integration
        """
        pass

    def getIncoherentExpectationValue(self, saveStringParams, expecName = None):
        """
            This method is the heart of the expectationvalue class.
            It expects the getRAWObservable function to be implemented.
        """
        expecObservable = np.zeros(self.prsr.interpTime.size) 
        TBFpops   = self.psFile.getTBFpopulations()
        rawObservable = self.getRAWObservable()
        # The following list contains the expectation value
        # calculated via the incoherent sum method for each IC.
        # If multiple runs have to be done it will be a list of
        # lists containing the expectation value for each run.  
        geomExpects = []
        ind = 0 
        if hasattr(self.prsr, "sampleSize"):
            for geom in np.arange(1, self.prsr.sampleSize + 1):
                if geom in self.prsr.dupList:
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
                        tmpCWD  = self.CWD + "/" + self.prsr.RNGdir + str(rng + 1)  
                        tmpCWD += "/" + self.prsr.geomDir + str(geom)
                        saveString = ""
                        for saveStringParam in saveStringParams: 
                            saveString += saveStringParam 
                        if not(expecName == None):
                            saveString += expecName[ind]
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
                    saveString += "_IC" + str(geom) + ".dat"
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
                    tmpCWD = self.CWD + "/" + self.prsr.geomDir + str(geom)
                    saveString  = "" 
                    for saveStringParam in saveStringParams:
                        saveString += saveStringParam
                    saveString += ".dat"
                    saveFile    = self.psFile.inputFileName(tmpCWD, saveString)
                    writeNPFile(2, saveFile, [self.prsr.interpTime, tmpExpect],
                                fmtStyle = "%8d %30.18e")

            saveFile = self.CWD + "/" 
            for saveStringParam in saveStringParams:
                saveFile += saveStringParam

            if self.prsr.AIMStype == "AIMS":
                nrSamples = self.prsr.sampleSize - len(self.prsr.dupList)
                saveFile += "_AIMS.dat"
            else:
                nrSamples  = (self.prsr.sampleSize - len(self.prsr.dupList)) 
                nrSsamples *= self.prsr.nrRNGs
                saveFile += "_" + self.prsr.AIMStype + "_" 
                saveFile += str(self.prsr.nrRNGs) + ".dat"
            expecObservable = expecObservable / nrSamples
            stdErrObservable = np.zeros(expecObservable.size)
            for i in np.arange(len(geomExpects)):
                stdErrObservable += (geomExpects[i] - expecObservable)**2
            stdErrObservable = np.sqrt(stdErrObservable / (nrSamples * 
                               (nrSamples + 1)))
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

    def interpTraj(self, traj, numParticles):
        time = []
        for snapshot in traj: 
            time.append(snapshot[0])
        interpCoords = []
        for j in np.arange(numParticles):
            for i in np.arange(3):
                currCoord = []
                for snapshot in traj: 
                    #print len(snapshot[1])
                    currCoord.append(snapshot[1][j][i])
                #print currCoord
                #print time, currCoord
                interpCurrCoord = np.interp(self.prsr.interpTime, 
                                            time, currCoord).tolist()
                interpCoords.append(interpCurrCoord) 

        newTraj = [] 
        for i in np.arange(self.prsr.interpTime.size):
            currTraj = []
            for j in np.arange(0,len(interpCoords)-2,3): 
                currTraj.append([interpCoords[j][i],
                                 interpCoords[j+1][i],
                                 interpCoords[j+2][i]])
            newTraj.append(currTraj)
                
        return newTraj

    def sampleGeometry(self, TBFcoord, widths, atmNames, rng):
        randCoord = []  
        for i, coord in enumerate(TBFcoord):
            tmpRandCoord = [atmNames[i]]
            sigma = np.sqrt(1/(4 * widths[i])) 
            #print coord
            for j, tmpCoord in enumerate(coord):
                if ((self.prsr.internalType == "X") or
                    (self.prsr.internalType == "Y")):
                    if j > 0:
                #if np.abs(tmpCoord) <= 1.e-12:
                        tmpRandCoord.append(0.0)    
                        continue
                tmpRandCoord.append(rng.normal(loc=tmpCoord,
                                               scale=sigma))
            randCoord.append(tmpRandCoord)
        return randCoord

    def getBox(self, rngInternal):
        for iBox, box in enumerate(self.prsr.boxes):
            if ((box[0] <= rngInternal) and 
                (rngInternal < box[1])):
                break
        else:
            if (rngInternal < box[0]):
                iBox = -1
            elif (rngInternal >= box[-1]):
                iBox = self.prsr.nrBoxes

        return iBox

    def calcGaussianDensity(self, RNGcoord, TBFcoord, widths):
        exponent = 0
        prefactor = 1
        for i, coord in enumerate(TBFcoord):
            alpha = widths[i]  

            for j, tmpCoord in enumerate(coord):
                if ((self.prsr.internalType == "X") or
                    (self.prsr.internalType == "Y")):
                    if j > 0:
                        break
                exponent -= 2 * alpha * (tmpCoord - RNGcoord[i][j+1])**2
                prefactor *= np.sqrt(2 * alpha / np.pi)
        gaussianDensity = prefactor * np.exp(exponent)
        return gaussianDensity

    def calcComplexGaussian(self, RNGcoord, TBFcoord, TBFmom, 
                            TBFphase, widths):
        exponent = 0
        prefactor = np.exp(1.0j * TBFphase)
        for i, coord in enumerate(TBFcoord):
            alpha = widths[i]  

            for j, tmpCoord in enumerate(coord):
                if ((self.prsr.internalType == "X") or
                    (self.prsr.internalType == "Y")):
                    if j > 0:
                        break
                exponent -= alpha * (RNGcoord[i][j+1] - tmpCoord)**2
                exponent += 1.0j * TBFmom[i][j] * (RNGcoord[i][j+1] - tmpCoord)
                prefactor *= np.sqrt(np.sqrt(2 * alpha / np.pi))
        complexGaussian = prefactor * np.exp(exponent)
        return complexGaussian

    def addDensity(self, t, currDensity, currDensity2, 
                   rngCoord, TBFcoord, TBFmom, TBFphase, 
                   TBFpop, TBFamp, TBFstt, widths):

        support = 0
        density = 0
        for iTBF, pop in enumerate(TBFpop):
            if pop[t] > 1.e-12:
                tmpGaussDensity = self.calcGaussianDensity(rngCoord,
                                                           TBFcoord[iTBF][t],
                                                           widths)
                support += pop[t] * tmpGaussDensity
                density += pop[t] * tmpGaussDensity
                for jTBF in np.arange(0,iTBF):
                    if ((TBFstt[iTBF] == TBFstt[jTBF]) and
                       (TBFpop[jTBF][t] > 1.e-12)):
                        chiI = self.calcComplexGaussian(rngCoord,
                                                        TBFcoord[iTBF][t],
                                                        TBFmom[iTBF][t],
                                                        TBFphase[iTBF][t],
                                                        widths)
                        chiJ = self.calcComplexGaussian(rngCoord,
                                                        TBFcoord[jTBF][t],
                                                        TBFmom[jTBF][t],
                                                        TBFphase[jTBF][t],
                                                        widths)
                        addTerm = np.conj(TBFamp[iTBF][t] * chiI) 
                        addTerm *= TBFamp[jTBF][t] * chiJ
                        density += 2 * np.real(addTerm)

        term = density/support
        #if self.prsr.interpTime[t] > 600:
        #print support
        newDensity = currDensity + density/support
        newDensity2 = currDensity2 + (density/support)**2
        return newDensity, newDensity2

    def importanceSampling(self, currCWD, t, TBFcoord, TBFmom, 
                           TBFphase, TBFpop, TBFamp, TBFstt,
                           widths, atmNames):
        if hasattr(self.prsr, "IRandSeed"):
            rng = np.random.RandomState(self.prsr.IRandSeed)
        else:
            rng = np.random.RandomState()

        meanOld = np.zeros(self.prsr.nrBoxes) 
        stderrOld = np.zeros(self.prsr.nrBoxes) 
        bOld = 0
        meanNew = np.zeros(self.prsr.nrBoxes) 
        stderrNew = np.zeros(self.prsr.nrBoxes) 
        bNew = 0
        nrSamples = 0
        nrFldTriesBlw = 0
        nrFldTriesAbv = 0
        failed = False
        density = np.zeros(self.prsr.nrBoxes)
        density2 = np.zeros(self.prsr.nrBoxes)
        #densitySamples = []
        #for iB in np.arange(self.prsr.nrBoxes): 
        #    densitySamples.append([])
        tes = 0 
        for iTBF, pop in enumerate(TBFpop):
            #print pop[t]
            tes += pop[t] 
        #print tes

        while True:
            meanOld = meanNew.copy()
            stderrOld = stderrNew.copy()
            bOld = bNew
            eta = rng.uniform(0,tes)
            cumProb = 0 
           

            for iTBF, pop in enumerate(TBFpop):
                if pop[t] < 1.e-12:
                    continue
                cumProb += pop[t] 
                if eta <= cumProb:
                    rngCoord = self.sampleGeometry(TBFcoord[iTBF][t],widths,
                                                   atmNames,rng) 
                    rngInternal = self.calcRAWObservable(rngCoord)
                    b = self.getBox(rngInternal)
                    if b == -1:
                        nrFldTriesBlw += 1
                        failed = True
                        break
                    if b == self.prsr.nrBoxes:
                        nrFldTriesAbv += 1
                        failed = True
                        break
                    density[b], density2[b] = self.addDensity(
                                              t, density[b],
                                              density2[b],
                                              rngCoord, TBFcoord,
                                              TBFmom, TBFphase,
                                              TBFpop, TBFamp,
                                              TBFstt, widths
                                              )
                    nrSamples += 1
                    bNew = b
                    break

            if not(failed):
                meanNew = density / nrSamples
                stderrNew = np.sqrt(density2 / nrSamples - meanNew**2)
                    
                if nrSamples > 1:
                    changeDens  = np.linalg.norm(meanNew - meanOld)
                    changeDens /= np.linalg.norm(meanOld)
                    if ((changeDens < self.prsr.densThresh) and
                        not(bNew == bOld)):
                        break
            else:
                failed = False

        boxsize = self.prsr.boxes[0][1] - self.prsr.boxes[0][0]
        summ = 0
        for i in np.arange(self.prsr.boxes.shape[0]):
            summ += meanNew[i] 

        #assert abs(summ - 1.0) <= 1.e-12
        print summ
        #print nrSamples
        finalDensity = []
        for iBox, box in enumerate(self.prsr.boxes):
            for limit in box:
                finalDensity.append(meanNew[iBox]) 

        finalDensity = np.array(finalDensity)
        #xaxis = self.prsr.boxes.flatten()
        #saveFile = 'redDens_' + str(t) + '.dat'
        #writeNPFile(2, saveFile, [xaxis, finalDensity], 
        #            fmtStyle = "%30.18e %30.18e") 
        print nrSamples + nrFldTriesBlw + nrFldTriesAbv
        return finalDensity

    def calcOverlap(self, xi, xj, pi, pj, gammaI, gammaJ, 
                    alpha):
        deltaXij = xi - xj 
        deltaPij = pi - pj 
        real = -0.5 * alpha * deltaXij**2  
        real -= 0.125 * (deltaPij**2/alpha)
        centroid = 0.5 * (xi + xj) 
        imag = (pi*xi - pj*xj) - centroid * deltaPij
        exponent = real + 1.j * imag
        overlap = np.exp(exponent)
        return overlap
             
    def calcExactDensity(self, currCWD, t, TBFcoord, TBFmom, 
                         TBFphase, TBFpop, TBFamp, TBFstt,
                         widths, atmNames, maxDens):

        xaxis = np.unique(self.prsr.boxes.flatten())
        exactDens = np.zeros(xaxis.size)
        for iX, x in enumerate(xaxis):
            fullDens = 0 
            nonOv = 0
            for iTBF in np.arange(len(TBFamp)):
                if self.prsr.internalType == "X":
                    tmpCoord = [TBFcoord[iTBF][t][0]]
                    tmpWidths = [widths[0]]
                    currCoord = [['QQ1', x, 0.0, 0.0]]
                elif self.prsr.internalType == "Y":
                    tmpCoord = [TBFcoord[iTBF][t][1]]
                    tmpWidths = [widths[1]]
                    currCoord = [['PP2', x, 0.0, 0.0]]
                gaussianDensity = self.calcGaussianDensity(
                                                 currCoord,
                                                 tmpCoord,
                                                 tmpWidths
                                                 )
                exactDens[iX] += TBFpop[iTBF][t] * gaussianDensity 
                fullDens += TBFpop[iTBF][t]
                nonOv += TBFpop[iTBF][t]
                for jTBF in np.arange(0,iTBF):
                    if TBFpop[iTBF][t] < 1.e-12:
                        continue
                    if TBFpop[jTBF][t] < 1.e-12:
                        continue
                    if not(TBFstt[iTBF] == TBFstt[jTBF]):
                        continue
                    #if self.prsr.interpTime[t] == 200: 
                    #    print iTBF, TBFpop[iTBF][t], jTBF, TBFpop[jTBF][t]
                    #print iTBF + 1, jTBF + 1, TBFpop[iTBF][t], TBFpop[jTBF][t]
                    widthX = widths[0]
                    widthY = widths[1]
                    overlapX = self.calcOverlap(
                              TBFcoord[iTBF][t][0][0],
                              TBFcoord[jTBF][t][0][0],
                              TBFmom[iTBF][t][0][0],
                              TBFmom[jTBF][t][0][0],
                              TBFphase[iTBF][t],
                              TBFphase[jTBF][t],
                              widthX
                              ) 
                    overlapY = self.calcOverlap(
                              TBFcoord[iTBF][t][1][0],
                              TBFcoord[jTBF][t][1][0],
                              TBFmom[iTBF][t][1][0],
                              TBFmom[jTBF][t][1][0],
                              TBFphase[iTBF][t],
                              TBFphase[jTBF][t],
                              widthY
                              ) 
                    if self.prsr.internalType == "X":
                        overlap = overlapY
                        chiI = self.calcComplexGaussian(
                               currCoord,
                               [TBFcoord[iTBF][t][0]],
                               [TBFmom[iTBF][t][0]],
                               TBFphase[iTBF][t],
                               [widthX]
                               )
                        chiJ = self.calcComplexGaussian(
                               currCoord,
                               [TBFcoord[jTBF][t][0]],
                               [TBFmom[jTBF][t][0]],
                               TBFphase[jTBF][t],
                               [widthX]
                               )
                        
                    elif self.prsr.internalType == "Y":
                        overlap = overlapX
                        chiI = self.calcComplexGaussian(
                               currCoord,
                               [TBFcoord[iTBF][t][1]],
                               [TBFmom[iTBF][t][1]],
                               TBFphase[iTBF][t],
                               [widthY]
                               )
                        chiJ = self.calcComplexGaussian(
                               currCoord,
                               [TBFcoord[jTBF][t][1]],
                               [TBFmom[jTBF][t][1]],
                               TBFphase[jTBF][t],
                               [widthY]
                               )
                    addTerm = 0
                    addTerm = np.conj(TBFamp[iTBF][t]) * TBFamp[jTBF][t] 
                    addTerm *= overlap * np.conj(chiI) * chiJ
                    addTerm = 2 * np.real(addTerm)
                    add2Term = 2 * np.real(np.conj(TBFamp[iTBF][t]) * TBFamp[jTBF][t]
                                            * overlapX * overlapY* np.exp(1.j*(TBFphase[jTBF][t] - TBFphase[iTBF][t])))
                    #if self.prsr.interpTime[t] == 200: 
                    #    print add2Term
                    #    print np.conj(TBFamp[iTBF][t]), TBFamp[jTBF][t]
                    #    print iTBF, jTBF
                    fullDens += add2Term  
                    #if addTerm > 1.e-12:
                    #    print iTBF, jTBF, addTerm
                    exactDens[iX] += addTerm
                    #print addTerm

            #print fullDens, nonOv, nonOv - 1
        intgrl = integrate.simps(exactDens, xaxis)
        print 'int', intgrl
        saveFile = 'exactRedDens_' + str(t) + '.dat'
        #exactDens = exactDens * maxDens / np.amax(exactDens)
        writeNPFile(2, saveFile, [xaxis, exactDens], 
                    fmtStyle = "%30.18e %30.18e") 
            
        return exactDens
                         
                        
    def calcCoherentExpectationValue(self, tmpCWD, saveStringParams, TBFpop, TBFamp):
        numParticles = self.psFile.findNumAtoms(tmpCWD)
        spawnTimes, childIDs, parentIDs, numSpawns = self.psFile.findNrSpawns(tmpCWD)
        if ((self.prsr.internalType == "X") or
            (self.prsr.internalType == "Y")):
            posErr, FGcoord = self.psFile.readPositions(1,tmpCWD,numParticles,
                                                        addAtmNames=False)
        else:
            posErr, FGcoord = self.psFile.readPositions(1,tmpCWD,numParticles,
                                                        addAtmNames=False, 
                                                        bohr=True)
        momErr, FGmom = self.psFile.readMomenta(1,tmpCWD,numParticles,
                                                addAtmNames=False)
        phaseErr, FGphase = self.psFile.getTBFphase(1,tmpCWD)
        TBFcoord = []
        TBFmom = []
        TBFphase = [] 
        if not(posErr or momErr or phaseErr):
            FGcoord = self.interpTraj(FGcoord, numParticles)
            FGmom = self.interpTraj(FGmom, numParticles)
            TBFcoord.append(FGcoord)
            TBFmom.append(FGmom)
            TBFphase.append(FGphase)
        for childID in childIDs:
            if ((self.prsr.internalType == "X") or
                (self.prsr.internalType == "Y")):
                posErr, CHcoord = self.psFile.readPositions(childID, 
                                                            tmpCWD,
                                                            numParticles,
                                                            addAtmNames=
                                                            False
                                                            )
            else:
                posErr, CHcoord = self.psFile.readPositions(childID, 
                                                            tmpCWD,
                                                            numParticles,
                                                            addAtmNames=
                                                            False,
                                                            bohr=True
                                                            )
            posErr, CHmom = self.psFile.readMomenta(childID, 
                                                    tmpCWD,
                                                    numParticles,
                                                    addAtmNames=
                                                    False)
            phaseErr, CHphase = self.psFile.getTBFphase(childID,
                                                        tmpCWD)
            if not(posErr or momErr or phaseErr):
                CHcoord = self.interpTraj(CHcoord, numParticles)
                CHmom = self.interpTraj(CHmom, numParticles)
                TBFcoord.append(CHcoord)
                TBFmom.append(CHmom)
                TBFphase.append(CHphase)
        TBFstt = self.psFile.getTBFstate(tmpCWD)
        widths = self.psFile.readWidths(tmpCWD) 
        atmNames = self.psFile.readAtomNames(tmpCWD) 
        redDens = []
        exactRedDens = []
        yaxis = self.prsr.boxes.flatten()
        densMovie = np.zeros((self.prsr.interpTime.size,yaxis.size,1))
        exactDensMovie = np.zeros((self.prsr.interpTime.size,yaxis.size,1))
        #print TBFmom
        for t in np.arange(self.prsr.interpTime.size):
            print self.prsr.interpTime[t]
            density = self.importanceSampling(tmpCWD, t, TBFcoord, TBFmom,
                                              TBFphase, TBFpop, TBFamp, 
                                              TBFstt, widths, atmNames)
            redDens.append(density)
            if ((self.prsr.internalType == "X") or
                (self.prsr.internalType == "Y")):
                exactDensity = self.calcExactDensity(tmpCWD, t, TBFcoord, 
                                                     TBFmom, TBFphase, 
                                                     TBFpop, TBFamp,
                                                     TBFstt, widths, 
                                                     atmNames, 
                                                     np.amax(density))
                exactRedDens.append(exactDensity)

        for t in np.arange(self.prsr.interpTime.size):
            for y in np.arange(yaxis.size):
                densMovie[t, y, 0] = redDens[t][y] 
                if ((self.prsr.internalType == "X") or
                    (self.prsr.internalType == "Y")):
                    exactDensMovie[t, y, 0] = exactRedDens[t][y]

        header = ['Time (atu)', 'internal', 'red. density']
        grid = [self.prsr.interpTime, yaxis]
        fileName = tmpCWD + '/redDens' 
        for saveString in saveStringParams: 
            fileName += saveString
        fileName += '.dat'
        writeGridFile(fileName, grid, densMovie, 1, header)

        if ((self.prsr.internalType == "X") or
            (self.prsr.internalType == "Y")):
            header = ['Time (atu)', 'internal', 'red. density']
            fileName = 'exactRedDens.dat' 
            writeGridFile(fileName, grid, exactDensMovie, 1, header)

        return densMovie
        
        #for iT, t in enumerate(self.prsr.interpTime):
        #     densMovie[0, iT*xaxis.size:(iT+1)*xaxis.size].size = t
        #     densMovie[1, iT*xaxis.size:(iT+1)*xaxis.size].size = xaxis 
        #     densMovie[2, iT*xaxis.size:(iT+1)*xaxis.size].size = density[iT]

    def getCoherentExpectationValue(self, saveStringParams):
        yaxis = self.prsr.boxes.flatten()
        fullDensity = np.zeros((self.prsr.interpTime.size,yaxis.size,1))
        nrSamples = 0
        TBFpop = self.psFile.getTBFpopulations()
        TBFamp = self.psFile.getTBFamplitude()
        if hasattr(self.prsr, "sampleSize"):
            ind = 0
            for geom in np.arange(1, self.prsr.sampleSize + 1):
                if geom in self.prsr.dupList:
                    continue
                if self.prsr.AIMStype != "AIMS":
                    geomExpecObservable = np.zeros(self.prsr.interpTime.size)
                    tmpExpects = []  
                    for rng in np.arange(self.prsr.nrRNGs):
                        print geom, rng
                        tmpCWD  = self.CWD + "/" + self.prsr.RNGdir 
                        tmpCWD += str(rng + 1) +  "/" + self.prsr.geomDir
                        tmpCWD += str(geom)
                        tmpDensity = self.calcCoherentExpectationValue(
                                        tmpCWD, saveStringParams, 
                                        TBFpop[ind][rng], TBFamp[ind][rng]
                                     )
                        fullDensity += tmpDensity  
                        nrSamples += 1
                    ind += 1
                else:
                    print geom
                    tmpCWD  = self.CWD + "/" + self.prsr.geomDir + str(geom)
                    tmpDensity = self.calcCoherentExpectationValue(
                                    tmpCWD, saveStringParams, 
                                    TBFpop[ind], TBFamp[ind]
                                 )
                    fullDensity += tmpDensity  
                    nrSamples += 1
                    ind += 1
        else:
            tmpCWD = self.CWD 
            fullDensity = self.calcCoherentExpectationValue(
                          tmpCWD, saveStringParams,
                          TBFpop, TBFamp
                      )
        
        fullDensityMean = fullDensity / nrSamples
        header = ['Time (atu)', 'internal', 'red. density']
        grid = [self.prsr.interpTime, yaxis]
        fileName = self.CWD + '/redDens' 
        for saveString in saveStringParams: 
            fileName += saveString
        fileName += '.dat'
        writeGridFile(fileName, grid, fullDensityMean, 1, header)

class internals(observables):
    def __init__(self, parser, cwd, psFile):
        observables.__init__(self, parser, cwd, psFile)
        self.observableType = "internals" 

    def getInternals(self):
        internals = 0
        if self.prsr.expecType == 'incoherent':
            self.getIncoherentExpectationValue(
                             saveStringParams = 
                             [self.prsr.internalType,
                             self.prsr.internalName]
                             )
        elif self.prsr.expecType == 'coherent':
            self.getCoherentExpectationValue(
                           saveStringParams = 
                           [self.prsr.internalType,
                           self.prsr.internalName]
                           )

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
        if (internalName == None):
            FGinternal, FGtime = self.calcInternals(FGcoord)
        else:
            FGinternal, FGtime = self.calcInternals(FGcoord, 
                                                    internalName =
                                                    internalName)
        interpFGintrnl = np.interp(self.prsr.interpTime, FGtime, 
                                   FGinternal)
        if (internalName == None):
            internal.append(interpFGintrnl)
        else:
            internalTime.append(FGtime[FGinternal > 0])
            internal.append(FGinternal[FGinternal > 0])

        if (internalName == None):
            self.saveInternals(tmpCWD, interpFGintrnl, 1)
         
        for childID in childIDs:
            posErr, CHcoord = self.psFile.readPositions(childID, 
                                                        tmpCWD,
                                                        numParticles)
            if not(posErr):
                if (internalName == None):
                    CHinternal, CHtime = self.calcInternals(CHcoord)
                else:
                    CHinternal, CHtime = self.calcInternals(CHcoord, 
                                                            internalName =
                                                            internalName)
                interpCHintrnl = np.interp(self.prsr.interpTime,  
                                           CHtime, CHinternal)
                if (internalName == None):
                        self.saveInternals(tmpCWD, interpCHintrnl, childID)

                if (internalName == None):
                    internal.append(interpCHintrnl)
                else:
                    internalTime.append(CHtime[CHinternal > 0])
                    internal.append(CHinternal[CHinternal > 0])

    def getRAWObservable(self): 
        internal = []
        if hasattr(self.prsr, "sampleSize"): 
            #if not(self.internalName == None):
            #    ind = 0
            for geom in np.arange(1, self.prsr.sampleSize + 1):
                geomInternal = []
                if geom in self.prsr.dupList:
                    continue
                if self.prsr.AIMStype != "AIMS":
                    for rng in np.arange(1, self.prsr.nrRNGs + 1):
                        rngInternal = []
                        tmpCWD = self.CWD + "/" + self.prsr.RNGdir + str(rng)  
                        tmpCWD += "/" + self.prsr.geomDir + str(geom)
                        self.addInternal(rngInternal, tmpCWD)
                        geomInternal.append(rngInternal)
                else:
                    tmpCWD = self.CWD + "/" + self.prsr.geomDir + str(geom)
                    self.addInternal(geomInternal, tmpCWD) 
                if len(geomInternal) != 0:
                    internal.append(geomInternal)
        else:
            tmpCWD = self.CWD 
            self.addInternal(internal, tmpCWD)

        return internal

    def calcRAWObservable(self, molStruct):
        atmsInvolved = self.prsr.internalName.split("-") 
        if self.prsr.internalType == "bl":
            assert (len(atmsInvolved) == 2)
        elif self.prsr.internalType == "ba":
            assert (len(atmsInvolved) == 3)
        elif self.prsr.internalType == "td":
            assert (len(atmsInvolved) == 4)
        elif ((self.prsr.internalType == "X") or
              (self.prsr.internalType == "Y")):
            assert (len(atmsInvolved) == 1)

        atmcoords = []
        for atm in molStruct:
            for i in np.arange(len(atmsInvolved)):
                if atm[0] == atmsInvolved[i]:
                    if not((self.prsr.internalType == "X") or
                           (self.prsr.internalType == "Y")):
                        atmcoords.append(np.array(atm[1:])*b2A)
                    else:
                        atmcoords.append(np.array(atm[1]))

        if not((self.prsr.internalType == "X") or
               (self.prsr.internalType == "Y")):
            rawObservable = calcIntrnl(atmcoords, self.prsr.internalType)
        else: 
            rawObservable = atmcoords[0]
        
        return rawObservable

