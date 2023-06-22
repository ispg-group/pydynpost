#!/usr/bin/env python
import numpy as np
import src.aims.inp as aimsinp
import src.tsh.inp as tshinp


pckgQuestions = {'aims': aimsinp.pckgQuestion,
                 'tsh': tshinp.pckgQuestion}

def initParser(parser):
    try:
        parser.addInput("pckg", pckgQuestions[parser.dynMethod])
    except:
        raise ValueError(f'Dynamics method {parser.dynMethod} not' +
                         ' recognized.')
    parser.addInput("todo", "What do you want to do?", arr = True)
    parser.addInput("geomDir", "What is the common IC directory name?")
    if not hasattr(parser, "geomDir"): 
        parser.addInternal("geomDir", "geom")
    if "population" in parser.todo:
        parser.addInput("nrStates", "How many states?")
    parser.addInput("model", "Model system?")
    if not hasattr(parser, "model"):
        parser.addInternal("model", parser.pckg)
    parser.addInput("nrStates", "How many states?")
    parser.addInput("nrParticles", "How many states?")
    parser.addInput("maxTime", "What is the largest time step?")
    parser.addInput("step", "Which resolution do you want?")
    parser.addInternal("interpTime", np.arange(0, parser.maxTime + parser.step, 
                                        parser.step))
    
    parser.addInput("sampleSize", "How many ICs?")
    parser.addInput("duplicatesExist", "Do duplicates exist?")
    if parser.duplicatesExist in ["y", "yes"]:
        parser.addInput("dupList", "Duplicates:", arr = True)
    else: 
        parser.addInternal("dupList", [])
    parser.addInput("multipleRNGs", "Are there multiple runs per IC" +
                                    " with different RNG seed?")
    if parser.multipleRNGs in ["y", "yes"]:
        parser.addInput("nrRNGs", "How many rngs?")
        parser.addInput("RNGdir", "What is the common RNG directory name?")
        if not hasattr(parser, "RNGdir"): 
            parser.addInternal("RNGdir", "rng")

    if "molpop" in parser.todo:
        parser.addInput("internalType", "Which kind of internal is it?")
        parser.addInput("dissPartners", "Which kind of internal is it?")
        parser.addInput("thresh", "Which kind of internal is it?")
        assert (hasattr(parser, "internalType") and hasattr(parser, "dissPartners")) 

    if "internals" in parser.todo:
        parser.addInput("internalType", "Which kind of internal is it?")
        parser.addInput("internalName", "Which atoms contribute?")
        assert (hasattr(parser, "internalType") and hasattr(parser, "internalName")) 

    parser.addInput("parallel","Parallel scripting?")
    parser = aimsinp.initParser(parser)

    return parser
