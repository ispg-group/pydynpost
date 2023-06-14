#!/usr/bin/env python
import numpy as np
import os
from commonmethods.filesys import *
from commonmethods.misc import *
from commonmethods.parse import *
from .aimsinp import *

def writeNPFile(nrColumns, fileName, content, fmtStyle = None):
    saveFormat      = np.zeros((content[0].size, nrColumns))
    for i in np.arange(nrColumns):
        saveFormat[:,i] = content[i]

    if fmtStyle != None:
        np.savetxt(fileName, saveFormat, fmt=fmtStyle)
    else:
        np.savetxt(fileName, saveFormat)

def writeGridFile(fileName, grid, scalars, numScalars, header,
                  width, precision):
    with open(fileName, 'w') as gridFile: 
        gridFile.write('# ')
        for columnDscrptn in header:
            fmt = ''.join(['{dscrptn:', width,'s}'])
            gridFile.write(fmt.format(dscrptn=columnDscrptn))
        gridFile.write('\n')

        for xi, x in enumerate(grid[0]):
            for yi, y in enumerate(grid[1]):
                fmt = '  ' + ''.join(['{x:', width, '.', precision, 'f} '])
                fmt += ''.join(['{y:', width, '.', precision,'f}']) 
                gridFile.write(fmt.format(x=x, y=y))
                for i in  np.arange(numScalars):
                    fmt = ''.join([' {scalar:', width, '.', precision, 'f}']) 
                    gridFile.write(fmt.format(scalar=scalars[xi, yi, i]))
                gridFile.write('\n')
            gridFile.write('\n')

               
