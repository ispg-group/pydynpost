#!/usr/bin/env python
import numpy as np

def writeNPFile(fileName, nrColumns, content, fmtStyle = None):
    saveFormat      = np.zeros((content[0].size, nrColumns))
    for i in np.arange(nrColumns):
        saveFormat[:,i] = content[i]

    if fmtStyle != None:
        np.savetxt(fileName, saveFormat, fmt=fmtStyle)
    else:
        np.savetxt(fileName, saveFormat)

def writeGridFile(fileName, grid, scalars, numScalars, header):
    with open(fileName, 'w') as gridFile: 
        gridFile.write('# ')
        for columnDscrptn in header:
            gridFile.write('{dscrptn:13s}'.format(dscrptn=
                                                  columnDscrptn))
        gridFile.write('\n')

        for xi, x in enumerate(grid[0]):
            for yi, y in enumerate(grid[1]):
                gridFile.write('  {x:13.6f} {y:13.6f}'.format(x=x, y=y))
                for i in  np.arange(numScalars):
                    gridFile.write(' {scalar:13.6f}'.format(scalar=
                                                            scalars[xi,
                                                                    yi,
                                                                    i]))
                gridFile.write('\n')
            gridFile.write('\n')
