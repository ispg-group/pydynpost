#!/usr/bin/env python
import os
import importlib
import numpy as np

class globalClass(object):
    def __init__(self, glbl):
        self.glbl = glbl

    def __getattr__(self, name):
        if hasattr(self.glbl, name):
            return getattr(self.glbl, name)
        elif name in self.__dict__:
            return self.__dict__[name]
        else:
            raise AttributeError(
                'Attribute not found in {0} and global class!'.format(
                    self.__class__.__name__
                )
            )
