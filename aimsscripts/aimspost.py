#!/usr/local/Cluster-Apps/python/2.7.9/bin/python

class AIMSpost(object):
    """
       Main postprocessing script for FMS/AIMS simulations. Gathers different 
       kinds of data and does statistics on it:
         - State Populations  
         - Number of TBFs
         - Theoretical Nr. of ES Calls (Worst Case)
         - Expectation values calculated in different ways
            *) Incoherent
            *) Coherent (not yet implemented, shipped at a later date maybe)
       Other miscellaneous methods include:
         - Calculation of momentum difference between parent and child at spawn 
         - Failure times 
         - Statistics of consistent histories (still to be finished)
       It also supports calculation of various decoherence times (ONLY for AIMS).
       Currently implemented:
         - Schwartz-Rossky Decoherence Time
    """
    moduleNames = modandmet.moduleNames
    methodNames = modandmet.methodNames

    def __init__(self, parser, cwd):
        # Important variables that almost all methods need to function
        self.prsr = parser
        self.CWD  = cwd 
        self.dirsInCwd = getDirs(self.CWD)
        self.psFile = processFiles(self.prsr, self.CWD)

        # Next only import those modules that are actually required
        # for a given task  
        importedModules = []
        for task in self.prsr.todo:
            moduleName = self.moduleNames[task] 
            # Do not load the same module twice 
            if moduleName in importedModules: 
                continue
            newModule = importlib.import_module(self.moduleNames[task])
            newClass  = getattr(newModule, self.methodNames[task])
            if task in ["population", "complexity"]: 
                # some tasks need to know the directory structure
                setattr(self, self.methodNames[task], 
                              newClass(self.prsr, self.CWD, 
                              self.dirsInCwd, self.psFile))
            else:
                setattr(self, self.methodNames[task], 
                              newClass(self.prsr, self.CWD, 
                              self.psFile))

    def performTasks(self):
        for task in self.prsr.todo:
            getterName = "get" + task[0].upper() + task[1:]
            currentTask = getattr(getattr(self, self.methodNames[task]), getterName)
            currentTask()


def main():
    cwd = os.getcwd()
    parser = parseInput("aimspost.inp", cwd)
    parser = initParser(parser) 
    postprocessing = AIMSpost(parser, cwd)
    postprocessing.performTasks()

if __name__ == "__main__":
    import numpy as np
    from matplotlib import pyplot as plt
    import os
    from matplotlib import cm
    import matplotlib as mpl
    import math
    import importlib
    import modandmet
    from filesys import *
    from misc import *
    from parse import *
    from aimsinp import *
    from aimsfiles import *
    main()
