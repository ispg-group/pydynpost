#!/usr/bin/env python
import numpy as np
import zarr


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

def getChunkShape(data):
    chunk = 10
    if data.shape[0] > data.shape[1]:
        while True:
            testarray = np.ones((chunk,data.shape[1]),dtype=data.dtype)
            if testarray.nbytes > int(1e6):
                break
            else:
                chunk *= 10

        return (chunk, data.shape[1])
    else:
        while True:
            testarray = np.ones((data.shape[0],chunk),dtype=data.dtype)
            print(testarray.size)
            if testarray.nbytes > int(1e6):
                break
            else:
                chunk *= 10

        return (data.shape[0], chunk)


def writeDataset(groupName, datasetName, content, parallel = False):
    chunkShape = getChunkShape(content) 
    dataset = groupName.create_dataset(
        datasetName, shape=content.shape, 
        chunks=chunkShape, dtype=content.dtype
    )
    if parallel == False:
        dataset[:] = zarr.array(
            content, shape=content.shape, 
            chunks=chunkShape, dtype=content.dtype
        )
    else:
        import dask.array
        groupName = groupName.__dict__['_path']
        dsName = "/" + groupName + "/" + datasetName
        content=content.rechunk(chunkShape)
        dask.array.to_zarr(content,'./data/griddata.zarr',
                           dsName, overwrite=True)

    return dataset.info 

