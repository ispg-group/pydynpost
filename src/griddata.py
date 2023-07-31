#!/usr/bin/env python
import importlib

def getGridData(glbl, geom, rng=None):
    readModule = 'src.' + glbl.dynMethod + '.files'
    readMethods = importlib.import_module(readModule)
    parseCheckpointFileFor = getattr(readMethods, 'parseCheckpointFileFor')
    err, outData = parseCheckpointFileFor(glbl, ['pgrid', 'pot', 'cgrid', 'coup'],
                                          geom, rng)
    outData = (outData['pgrid'], outData['pot'], outData['cgrid'], outData['coup']) 
    if err == True:
        raise FileNotFoundError(f"Checkpoint file for IC ${geom} not found.")
    else:
        return outData
