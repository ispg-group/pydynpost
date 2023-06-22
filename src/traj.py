#!/usr/bin/env python
import os
import importlib
import itertools
import numpy as np
import globalattr as glAttr

class Traj(glAttr.globalClass):

    def __init__(self, glbl, geom, trajID, rng = None):

        super().__init__(glbl)
        self.geom = geom
        self.trajID = trajID
        self.rng = rng

    def getPositions(self):
        readMethods = importlib.import_module(dynMethod + 'files')
        readPositions = getattr(readMethods, 'readPositions')
        err, positions = readPositions(
            self, self.trajID, self.geom, rng = self.rng
        )
        if not(err):
            self.positions = postions
        else:
            raise FileNotFoundError(err)

    def getMomenta(self):
        readMethods = importlib.import_module(dynMethod + 'files')
        readMomenta = getattr(readMethods, 'readMomenta')
        err, momenta = readMomenta(
            self, self.trajID, self.geom, rng = self.rng
        )
        if not(err):
            self.momenta = momenta
        else:
            raise FileNotFoundError(err)


    def getPhase(self):
        if self.dynMethod == 'tsh':
            self.calculatePhase()
        else:
            readMethods = importlib.import_module('aimsfiles')
            readPhase = getattr(readMethods, 'readPhase')
            err, phase = trajModule.getPhase(
                self, self.trajID, self.geom, rng = self.rng
            )
            if not(err):
                self.momenta = momenta
            else:
                raise FileNotFoundError(err)
            

    def calculatePhase(self):
        raise NotImplementedError

    def getAmplitude(self):
        if self.dynMethod == 'tsh':
            return True 
        else:
            readMethods = importlib.import_module('aimsfiles')
            readAmplitude = getattr(readMethods, 'readAmplitude')
            err, amplitude = readAmplitude(
                self, self.trajID, self.geom, rng = self.rng
            )

            if not(err):
                self.amp = amplitude
            
            return err

    def getAmplitudeDot(self):
        if self.dynMethod == 'tsh':
            return True 
        else:
            readMethods = importlib.import_module('aimsfiles')
            readAmplitudeDot = getattr(readMethods, 'readAmplitudeDot')
            err, amplitudeDot = readAmplitudeDot(
                self, self.trajID, self.geom, rng = self.rng
            )

            if not(err):
                self.ampDot = amplitudeDot
            
            return err

def getSimpleIterator(glbl, geom, rng = None):

    nTraj = readNrTBFs(geom, rng = rng)[-1][-1]
    def trajSimpleIterator():

        for iTraj in range(1, nTraj + 1):

            currTraj = Traj(glbl, geom, itraj, rng = rng)
            yield currTraj
            
    simpleIterator = trajSimpleIterator()
         
    return simpleIterator 
