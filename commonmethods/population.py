#!/usr/bin/env python
import importlib

def getPopulation(dynMethod, pkcg, CWD, geom, geomDir, 
                  rng=None, rngDir=None):
    readMethods = importlib.import_module(dynMethod + 'files')
    readPopulation = getattr(readMethods, 'readPopulation')

    if not((rng == None) or (rngDir == None)):
        population = readPopulation(pkcg, CWD, geom, geomDir, rng, rngDir) 
    else:
        population = readPopulation(pkcg, CWD, geom, geomDir) 

    if dynMethod == 'tsh':
        import tshpopulation
        population = tshpopulation.getPopulation(population) 

    return population
