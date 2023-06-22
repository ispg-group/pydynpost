#!/usr/bin/env python

pckgQuestion = "Which AIMS interface did you use?"

def initParser(parser):
    if parser.pckg == "molpro":
        parser.addInput("outputDir", "What is the directory name?")

    return parser
