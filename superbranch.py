#!/usr/bin/env python
import os
import importlib
import numpy as np
import globalattr as glAttr
import branch

class superbranch(glAttr.globalClass):
    def __init__(self, glbl):
        super().__init__(glbl)
        self.branches = branch.getSimpleIterator(glbl)

    def getStatisticalTimeTrace(self, metric):
        partialSum = 0.0
        partialSquareSum = 0.0

        for branch in self.branches:
            tmpSummand = branch.getMetric(metric)
        #    tmpSquareSummand = tmpSummand ** 2 
        #    partialSum += tmpSummand
        #    partialSquareSum += tmpSquareSummand

        #if hasattr(self, 'nrRNGs'):
        #    totNrSamples = self.nrUnqSamples * self.nrRNGs 
        #else:
        #    totNrSamples = self.nrUnqSamples

        #mTimeTrace = partialSum / totNrSamples 
        #stdTimeTrace = partialSquareSum / totNrSamples - mTimeTrace**2 
        #stdTimeTrace = np.sqrt(1./(totNrSamples - 1) * stdTimeTrace)

        #return mTimeTrace, stdTimeTrace
    
    def propagatePost(self):
        pass

if __name__ == '__main__':
    import parse
    cwd = os.getcwd()
    parser = parse.parseInput("dynpost.inp", cwd)
    parser.addInput("dynMethod", "Which dynamics method was used?")
    inpName = parser.dynMethod + "inp"
    inpModule = importlib.import_module(inpName)
    parser = inpModule.initParser(parser) 
    pcFMS = superbranch(parser)
    pcFMS.getStatisticalTimeTrace('population')
