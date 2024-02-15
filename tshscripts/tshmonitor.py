#!/usr/bin/env python
import numpy as np
import os
import math
from commonmethods.filesys import *
from commonmethods.misc import *
from commonmethods.parse import *

class Monitor(object):
    def __init__(self, parser, cwd, psFile):
        self.prsr = parser
        self.CWD  = cwd
        self.psFile = psFile

    def getMonitor(self):
