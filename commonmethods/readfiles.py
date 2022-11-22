#!/usr/bin/env python
import numpy as np

def readNPfile(fileName, keyIndices):
    readData = np.genfromtxt(fileName) 
    outDict = {}
    for key, iKey in keyIndices.items(): 
        outDict[key] = readData[:, iKey] 

def readTraj(fileName, nrAtoms, readTimestep):
    """
        Since most trajectory files have the same
        structure but a different position of the 
        timestep, this must be supplied by the 
        calling procedure as the lambda function 
        readTimestep.
    """
    with open(fileName, "r") as pLines:
        readNrParticles = True
        startRead = False
        nrLines = 0
        posTimes = []
        positions = []
        curPos = []
        for pLine in pLines:
            curSpltLine = pLine.strip().split()
            if (len(curSpltLine) == 1) and (readNrParticles == True):
                nrParticles     = int(curSpltLine[0])
                readNrParticles = False 
                startRead = True
                nrLines += 1
                continue
            elif (len(curSpltLine) == 1): 
                startRead = True
                nrLines += 1
                continue

            if (startRead and (nrLines >= 2)):
                if (nrLines < (nrParticles + 1)):
                    atomName = curSpltLine[0]
                    atomPos  = [atomName + str(nrLines - 1)]
                    atomPos.extend([float(x) for x in curSpltLine[1:]])
                    curPos.append(atomPos)
                    nrLines += 1
                else:
                    atomName = curSpltLine[0]
                    atomPos  = [atomName + str(nrLines - 1)]
                    atomPos.extend([float(x) for x in curSpltLine[1:]])
                    curPos.append(atomPos)
                    positions.append(curPos)
                    startRead = False
                    nrLines = 0
                    curPos = []
            else:
                nrLines += 1
                curTime = readTimestep(curSpltLine) 
                posTimes.append(curTime)
                continue

    return list(zip(posTimes,positions))
