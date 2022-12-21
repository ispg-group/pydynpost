#!/usr/bin/env python
import os
import importlib
import numpy as np
import globalattr as glAttr
import branch

class Tree(glAttr.globalClass):

    def __init__(self, glbl):

        super().__init__(glbl)
        self.branches = branch.getSimpleIterator(glbl)

    def getIBExpectation(self, metric):

        partialSum = 0.0
        partialSquareSum = 0.0

        for branch in self.branches:

            summand = branch.getMetric(metric)
            squareSummand = summand ** 2 

            for leaf in branch.leaves:
                tmpSummand = leaf.getMetric(metric) 
                summand += tmpSummand
                squareSummand += tmpSummand ** 2 


            partialSum += summand
            partialSquareSum += squareSummand

        if hasattr(self, 'nrRNGs'):
            totNrSamples = self.nrUnqSamples * self.nrRNGs 

        totNrSamples = self.nrUnqSamples

        mTimeTrace = partialSum / totNrSamples 
        stdTimeTrace = partialSquareSum / totNrSamples - mTimeTrace**2 
        stdTimeTrace = np.sqrt(1./(totNrSamples - 1) * stdTimeTrace)

        return mTimeTrace, stdTimeTrace
    
    def propagatePost(self):

        pass

if __name__ == '__main__':
    import parse
    cwd = os.getcwd()
    parser = parse.parseInput("dynpost.inp", cwd)
    parser.addInput("dynMethod", "Which dynamics method was used?")
    import commoninp
    parser = commoninp.initParser(parser) 
    pcFMS = Tree(parser)
    pcFMS.getIBExpectation('population')
