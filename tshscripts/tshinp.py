#!/usr/bin/env python
import numpy as np
import os
from commonmethods.filesys import *
from commonmethods.misc import *
from commonmethods.parse import *

def initParser(parser):
    parser.addInput("code", "Which tsh code did you use")
    parser.addInput("todo", "What do you want to do?", arr = True)

    parser.addInput("geomDir", "What is the common IC directrory name?")
    if not(hasattr(parser, "geomDir")): 
        parser.addd("geomDir", "geom")
    if (("population" in parser.todo) or ("coupling" in parser.todo)):
        parser.addInput("nrStates", "How many states?")
    parser.addInput("maxTime", "What is the largest time step?")
    parser.addInput("step", "Which resolution do you want?")
    parser.addd("interpTime", np.arange(0, parser.maxTime + parser.step, 
                                        parser.step))

    parser.addInput("nrRNGs", "How many rngs?")
    parser.addInput("RNGdir", "What is the common RNG directory name?")
    if not(hasattr(parser, "RNGdir")): 
        parser.addd("RNGdir", "rng")

    parser.addInput("sampleSize", "How many ICs?")
    parser.addInput("duplicatesExist", "Do duplicates exist?")
    if parser.duplicatesExist in ["y", "yes"]:
        parser.addInput("dupList", "Duplicates:", arr = True)
    else: 
        parser.addd("dupList", [])

    if "internals" in parser.todo:
        parser.addInput("internalType", "Which kind of internal is it?")
        parser.addInput("internalName", "Which atoms contribute?")
        assert (hasattr(parser, "internalType") and hasattr(parser, "internalName")) 
    elif "molpop" in parser.todo:
        parser.addInput("internalType", "Which kind of internal is it?")
        parser.addInput("dissPartners", "Which kind of internal is it?")
        parser.addInput("thresh", "test")
        assert (hasattr(parser, "internalType") and hasattr(parser, "dissPartners")) 
    if "coupling" in parser.todo:
        parser.addInput("couplingType", "What type of effective nac was used?")
        parser.addInput("CSThresh", "What type of effective nac was used?")

    if "monitor" in parser.todo: 
        parser.addInput("monitorObservable", "Which observable do you want to monitor?")

    return parser
