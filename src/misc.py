#!/usr/bin/env python
import numpy as np
import os

def isint(x):
    merr = False 
    if type(x) == int:
        merr = True
    elif type(x) == float:
        pass
    elif type(x) == bool:
        pass
    elif type(x) == str:
        try: 
            x = int(x)
            merr = True
        except:
            pass

    return merr

def isfloat(x):
    merr = False
    if type(x) == int:
        pass
    elif type(x) == float:
        merr = True
    elif type(x) == bool:
        pass
    elif type(x) == str:
        try: 
            x = float(x)
            merr = True 
        except:
            pass

    return merr

def isbool(x):
    merr = False
    if type(x) == int:
        pass
    elif type(x) == float:
        pass
    elif type(x) == bool:
        merr = True
    elif type(x) == str:
        if x in ["True", "False"]:
            merr = True

    return merr

def convtobool(x):
    if type(x) == bool:
        return x
    elif type(x) == str:
        if x == "True":
            return True
        elif x == "False":
            return False
        else:
            raise ValueError("Cannot convert arbitrary string to boolean")
    else:
        return bool(x)

def isstring(x):
    merr = False
    if type(x) == int:
        pass
    elif type(x) == float:
        pass
    elif type(x) == bool:
        pass
    elif type(x) == str:
        merr = True

    return merr

def islist(x):
    merr = False
    if type(x) == int:
        pass
    elif type(x) == float:
        pass
    elif type(x) == bool:
        pass
    elif type(x) == str:
        pass
    elif type(x) == list:
        merr = True

    return merr
