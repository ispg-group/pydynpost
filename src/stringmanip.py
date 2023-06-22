#!/usr/bin/env python
import numpy as np
import os
from src.misc import *

def splitline(line, splittype):                                                                                                                                                                                                                                       
    lineContents = line.split(splittype)[1]
    if len(line.split()) <= 3:
        if isint(lineContents):
            return int(line.split(splittype)[1].strip())
        elif isfloat(lineContents):
            return float(line.split(splittype)[1].strip())
        elif isbool(lineContents):
            return bool(line.split(splittype)[1].strip())
        elif isstring(lineContents):
            return line.split(splittype)[1].strip()
    else: 
        outptLine = []
        lineContents = line.split(splittype)[1]
        lineContents = lineContents.split(",")
        for lineContent in lineContents:
            if isint(lineContent):
                outptLine.append(int(lineContent))
            elif isfloat(lineContent):
                outptLine.append(float(lineContent))
            elif isbool(lineContent):
                outptLine.append(bool(lineContent))
            elif isstring(lineContent):
                outptLine.append(str(lineContent))

        return outptLine

