#!/usr/bin/env python

class const(object):
    def __setattr__(self, name, value):                                                                                                                                                                                                                             
        if name in self.__dict__:
            raise ValueError("cannot change a const attribute")
        self.__dict__[name] = value 
                                    
    def __delattr__(self, name):
        if name in self.__dict__:
            raise ValueError("cannot delete a const attribute")
        raise AttributeError("'{0}' object has no attribute '{1}'".format(
                             self.__class__.__name__, name))
