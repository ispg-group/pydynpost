#!/usr/bin/env python
import numpy as np
import os
import math
import importlib
from commonmethods.filesys import *
from commonmethods.misc import *
from commonmethods.parse import *
from .tshinp import *
from . import tshfiles 
from . import tshmodandmet

class TSHpost(object):
    """
       Main postprocessing script for TSH simulations. Gathers different 
       kinds of data and does statistics on it:
         - State Populations 
         - Expectation values calculated via incoherent sum over 
           independent trajectories
       More to come later on
    """
    moduleNames = tshmodandmet.moduleNames
    methodNames = tshmodandmet.methodNames

    def __init__(self, parser, cwd):
        # Important variables that almost all methods need to function
        self.prsr = parser
        self.CWD  = cwd 
        self.dirsInCwd = getDirs(self.CWD)
        # find processFiles class correspoding to correct TSH code
        psFileClass = getattr(tshfiles, "processFiles" + self.prsr.code.upper()) 
        self.psFile = psFileClass(self.prsr, self.CWD)

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
    main()
