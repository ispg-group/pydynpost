#!/usr/bin/env python
from __future__ import print_function 
import numpy as np
import os
import sys
import const
from stringmanip import * 
from misc import *

importedPromptTK = True
try:
    from prompt_toolkit import print_formatter_text as print
except:
    importedPromptTK = False 


class parseInput(object):
    def __init__(self, inputFileName, CWD):
        self.usrInput = const.const() 
        self.inpFName = inputFileName
        self.CWD      = CWD
        self.addInternal("CWD", CWD)

    def _formatParseQ(self, question, totL = 5):
        questionLen = len(question) 
        numtab = totL - questionLen / 8;
        for tab in np.arange(numtab):
            question += "\t"
        return question 

    def addInternal(self, name, value):                                                                                                                                                                                                                             
        self.usrInput.__setattr__(name, value)

    def addInput(self, variableName, question, arr = None, totL = 5):
        if self.inpFName in os.listdir(self.CWD):
            tmpFName = self.CWD + "/" + self.inpFName
            with open(tmpFName, "r") as tmpFLines:
                for line in tmpFLines:
                    if variableName in line:
                        sep = line.split()[1]
                        values = splitline(line, sep) 
                        if (arr == None) or islist(values):
                            self.usrInput.__setattr__(variableName, values) 
                        elif (arr != None): 
                            valuescp = []
                            valuescp.append(values)
                            self.usrInput.__setattr__(variableName, valuescp) 
        else:
            try:
            # Python2.x way
                values = raw_input(self._formatParseQ(question, totL = totL)).split(",")
            except:
            # Python3.x way
                values = input(self._formatParseQ(question, totL = totL)).split(",")
            if (len(values) == 1) and (arr == None):
                value = values[0]
                if isint(value):
                    self.usrInput.__setattr__(variableName, int(value))
                elif isfloat(value):
                    self.usrInput.__setattr__(variableName, float(value))
                elif isbool(value): 
                    self.usrInput.__setattr__(variableName, convtobool(value))
                elif isstring(value):
                    self.usrInput.__setattr__(variableName, value)
            else:
                valuescp = []
                for value in values:
                    if isint(value):
                        valuescp.append(int(value)) 
                    elif isfloat(value):
                        valuescp.append(float(value))
                    elif isbool(value): 
                        valuescp.append(coconvnvtobool(value))
                    elif isstring(value):
                        valuescp.append(value)
                self.usrInput.__setattr__(variableName, valuescp)

    def __getattr__(self, name):
        try: 
            return self.usrInput.__dict__[name]
        except:
            raise AttributeError("'{0}' object has no attribute '{1}'".format(
                                 self.__class__.__name__, name))
#        if hasattr(self.usrInput, name):
#            return self.usrInput.__dict__[name]
#        else:
#            raise AttributeError("'{0}' object has no attribute '{1}'".format(
#                                 self.__class__.__name__, name))
