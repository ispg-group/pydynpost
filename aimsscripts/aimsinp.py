#!/usr/bin/env python
import numpy as np
import os
from commonmethods.filesys import *
from commonmethods.misc import *
from commonmethods.parse import *

def initParser(parser):
    parser.addInput("pckg", "Which AIMS interface was used?")
    parser.addInput("todo", "What do you want to do?", arr = True)
    if parser.pckg == "molpro": 
        parser.addInput("outputDir", "What is the directory name?")

    parser.addInput("geomDir", "What is the common IC directrory name?")
    if not(hasattr(parser, "geomDir")): 
        parser.addd("geomDir", "geom")
    parser.addInput("AIMStype", "What AIMS type was used?")
    if (("population" in parser.todo) or ("coupling" in parser.todo)):
        parser.addInput("nrStates", "How many states?")
    parser.addInput("maxTime", "What is the largest time step?")
    parser.addInput("step", "Which resolution do you want?")
    parser.addd("interpTime", np.arange(0, parser.maxTime + parser.step, 
                                        parser.step))

    if parser.AIMStype != "AIMS":
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
        parser.addInput("expecType", "Incoherent or coherent expectation value?")
        if parser.expecType == "coherent":
            parser.addInput("expecMin", "Minimum expectation value?")
            parser.addInput("expecMax", "Maximum expectation value?")
            parser.addInput("nrBoxes", "Number of boxes?")
            interval = np.linspace(parser.expecMin,parser.expecMax,
                                   parser.nrBoxes+1)
            boxes = []
            for i in np.arange(parser.nrBoxes):
                boxes.append((interval[i], interval[i+1])) 
            parser.boxes = np.array(boxes)
            parser.addInput("densThresh", "Threshold for MC sampling?")

        assert (hasattr(parser, "internalType") and hasattr(parser, "internalName")) 
    elif "molpop" in parser.todo:
        parser.addInput("internalType", "Which kind of internal is it?")
        parser.addInput("dissPartners", "Which kind of internal is it?")
        parser.addInput("thresh", "Which kind of internal is it?")
        assert (hasattr(parser, "internalType") and hasattr(parser, "dissPartners")) 
    if "coupling" in parser.todo:
        parser.addInput("couplingType", "What type of effective nac was used?")
        parser.addInput("CSThresh", "What type of effective nac was used?")

    return parser
