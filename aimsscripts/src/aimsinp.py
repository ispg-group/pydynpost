#!/usr/bin/env python

pckgQuestion = "Which AIMS interface did you use?"

def initParser(parser):
    if "internals" in parser.todo:
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

    if "coupling" in parser.todo:
        if not hasattr(parser, "nrStates")):
            parser.addInput("nrStates", "How many states?")
        parser.addInput("couplingType", "What type of effective nac was used?")
        parser.addInput("CSThresh", "What type of effective nac was used?")

    return parser
