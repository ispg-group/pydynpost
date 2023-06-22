#!/usr/bin/env python
import os
import importlib
import numpy as np

class globalClass(object):
    def __init__(self, glbl):
        self.glbl = glbl

    def __getattr__(self, name):
        if name in self.__dict__:
            return self.__dict__[name]
        try: 
            return getattr(self.glbl, name)
        except:
            raise AttributeError(
                'Attribute {0} not found in {1} and global class!'.format(
                    name, self.__class__.__name__
                )
            )
        #if hasattr(self.glbl, name):
        #    return getattr(self.glbl, name)
        #elif name in self.__dict__:
        #    return self.__dict__[name]
        #else:
        #    raise AttributeError(
        #        'Attribute {0} not found in {1} and global class!'.format(
        #            name, self.__class__.__name__
        #        )
        #    )
