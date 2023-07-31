#!/usr/bin/env python
import importlib
import numpy as np
from scipy import interpolate

def getPopulation(glbl, geom, rng=None):
    readModule = 'src.' + glbl.dynMethod + '.files'
    readMethods = importlib.import_module(readModule)
    readPopulation = getattr(readMethods, 'readPopulation')

    if not (rng == None):
        time, population = readPopulation(glbl, geom, rng) 
    else:
        time, population = readPopulation(glbl, geom) 

    if glbl.dynMethod == 'tsh':
        import tshpopulation
        population = tshpopulation.getPopulation(population) 
    
    population = _interpolatePopulation(glbl.interpTime, population, time)

    return population

def _interpolatePopulation(interpTime, population, time):
    
    interpPopulation = np.zeros((interpTime.size,population.shape[1]))
    for state, statePop in enumerate(population.T):
        popDuring = interpolate.interp1d(
            time, statePop, bounds_error=False, 
            fill_value=(statePop[0], statePop[-1])
        ) 
        interpPopulation[:, state] = popDuring(interpTime) 
    
    return interpPopulation
