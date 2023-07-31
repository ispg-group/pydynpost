#!/usr/bin/env python
import os
import importlib
import itertools
import numpy as np
import src.globalattr as glAttr
from src.aims.files import readNrTBFs

class Traj(glAttr.globalClass):

    def __init__(self, glbl, geom, trajID = None, rng = None, time = None):

        super().__init__(glbl)
        self.geom = geom
        self.trajID = trajID
        self.rng = rng
        self.filesModule = importlib.import_module(
            'src.' + self.dynMethod + '.files'
        )
        self.time = time
        print(self.geom, self.time, self.trajID)

    def getStateID(self):
        readStateID = getattr(self.filesModule, 'readStateID')
        err, stateID = readStateID(
            self, self.geom, rng = self.rng, 
            trajID = self.trajID  
        )
        if not(err):
            self.stateID = stateID
        else:
            raise IOError(err)
        

    def getPositions(self):
        readPositions = getattr(self.filesModule, 'readPositions')
        err, positions = readPositions(
            self, self.geom, rng = self.rng, 
            trajID = self.trajID, time = self.time  
        )
        if not(err):
            self.positions = positions
            print(positions)
        else:
            raise IOError(err)

    def getMomenta(self):
        readMomenta = getattr(self.filesModule, 'readMomenta')
        err, momenta = readMomenta(
            self, self.geom, rng = self.rng, 
            trajID = self.trajID, time = self.time
        )
        if not(err):
            self.momenta = momenta
            print(self.momenta)
        else:
            raise IOError(err)

    def getPhase(self):
        if self.dynMethod == 'tsh':
            self.calculatePhase()
        else:
            readPhase = getattr(self.filesModule, 'readPhase')
            err, phase = readPhase(
                self, self.geom, rng = self.rng, 
                trajID = self.trajID, time = self.time
            )
            if not(err):
                self.phase = phase
                print(self.phase)
            else:
                raise IOError(err)

    def getWidths(self):
        if self.dynMethod == 'tsh':
            raise NotImplementedError
        else:
            readWidths = getattr(self.filesModule, 'readWidths')
            err, widths = readWidths(
                self, self.geom, rng = self.rng
            )
            if not(err):
                self.widths = widths
                print(self.widths)
            else:
                raise IOError(err)
            
    def calculatePhase(self):
        raise NotImplementedError

    def getAmplitude(self):
        if self.dynMethod == 'tsh':
            return True 
        else:
            readAmplitude = getattr(self.filesModule, 'readAmplitude')
            err, amplitude = readAmplitude(
                self, self.geom, rng = self.rng, 
                trajID = self.trajID, time = self.time
            )

            if not(err):
                self.amp = amplitude
                print(self.amp)
            else:
                raise IOError(err) 

    def getAmplitudeDot(self):
        if self.dynMethod == 'tsh':
            return True 
        else:
            readAmplitudeDot = getattr(self.filesModule, 'readAmplitudeDot')
            err, amplitudeDot = readAmplitudeDot(
                self, self.geom, rng = self.rng, 
                trajID = self.trajID, time = self.time
            )

            if not(err):
                self.ampDot = amplitudeDot
            
            return err

    def calcAmplitudeDot(self): 
        readAmplitude = getattr(self.filesModule, 'readAmplitude')
        ind = np.argwhere(self.interpTime==self.time)[0,0]
        time = np.array(
            [self.interpTime[ind-2], self.interpTime[ind-1],
             self.interpTime[ind+1], self.interpTime[ind+2]]  
        )
        reverse = False
        print(self.time)
        if self.time == self.interpTime[0]: 
            times = self.interpTime[:3] 
            weights = np.array([-1.5, 2.0, -0.5])  
            print(self.time)
        elif self.time == self.interpTime[1]: 
            times = np.array(
                [self.interpTime[0], self.interpTime[2]]
            )
            weights = np.array([-0.5, 0.5])  
        elif self.time == self.interpTime[-1]: 
            lastInd = self.interpTime.size
            slicer = slice(lastInd, lastInd-4, -1) 
            times = self.interpTime[slicer] 
            weights = np.array([1.5, -2.0, 0.5])  
            reverse = True
        elif self.time == self.interpTime[-2]: 
            times = np.array(
                [self.interpTime[-1], self.interpTime[-3]]
            )
            weights = np.array([-0.5, 0.5])  
            reverse = True

        err, amplitudes = readAmplitude(
            self, self.geom, rng = self.rng, 
            trajID = self.trajID, time = times,
            reverse = reverse
        )
        if err:
            raise IOError(err) 
        self.ampDot = np.dot(weights, amplitudes)  
        print(amplitudes, weights, self.ampDot)
    
    def pruneData(self): 
        for key in self.positions.keys():
            self.positions[key] = np.array([self.positions[key][0]])
            self.momenta[key] = np.array([self.momenta[key][0]])

        self.widths = [ self.widths[i] 
                        for i in range(0,3*self.nrParticles,3) ] 

class TrajSimpleIterator():
    def __init__(self, glbl, nTraj, geom, rng = None, time = None):
        self.glbl = glbl
        self.geom = geom
        self.rng = rng
        self.time = time
        self.nTraj = nTraj 

    def __iter__(self):
        self.iTraj = 0
        return self 

    def __next__(self):
        if self.iTraj < (self.nTraj):
            self.iTraj += 1
            return Traj(
                self.glbl, self.geom, trajID=self.iTraj, 
                rng=self.rng, time=self.time
            )

        else:
            raise StopIteration

def getSimpleIterator(glbl, geom, rng = None, time = None):
    if time != glbl.interpTime[0]:
        nTraj = readNrTBFs(glbl, geom, rng = rng)[-1][-1]
    else:
        nTraj = 1
        
    simpleIterator = TrajSimpleIterator(
        glbl, nTraj, geom, rng=rng, time=time
    )
         
    return nTraj, iter(simpleIterator)

def constructPairIterator(trajs1, trajs2 = None, diagonals = False):
    iTraj1 = 0 
    for traj1 in trajs1: 
        if trajs2 == None:
            jTraj1 = 0
            for traj in trajs1:
                if ((diagonals == True) and 
                    traj.trajID == traj1.trajID):
                    jTraj1 += 1 
                    continue
                yield (iTraj1, jTraj1), (traj1, traj) 
                jTraj1 += 1 
        else:
            jTraj2 = 0
            for traj2 in trajs2:
                print((iTraj1, jTraj2), (traj1, traj2))
                yield (iTraj1, jTraj2), (traj1, traj2) 
                jTraj2 += 1 

        iTraj1 += 1
