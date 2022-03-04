#!/usr/local/Cluster-Apps/python/2.7.9/bin/python
import numpy as np
from matplotlib import pyplot as plt
import os
from matplotlib import cm
import matplotlib as mpl
import math
from filesys import *
from misc import *
from parse import *
from aimsinp import *

def writeNPFile(nrColumns, fileName, content, fmtStyle = None):
    saveFormat      = np.zeros((content[0].size, nrColumns))
    for i in np.arange(nrColumns):
        saveFormat[:,i] = content[i]

    if fmtStyle != None:
        np.savetxt(fileName, saveFormat, fmt=fmtStyle)
    else:
        np.savetxt(fileName, saveFormat)
