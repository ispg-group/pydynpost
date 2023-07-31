#!/usr/bin/env python
import numpy as np
from matplotlib import pyplot as plt
import os
from matplotlib import cm
import matplotlib as mpl
import math

grad2rad = np.pi / 180.
rad2grad = 180. / np.pi
b2A = 0.529177249
A2b = 1./b2A

def calcBL(atmcoordi, atmcoordj):
    """
        Calculate the length between two vectors.
        Mostly used to calculate bond-length
    """
    #print atmcoordi, atmcoordj
    vecji = atmcoordi - atmcoordj
    return np.sqrt(np.dot(vecji, vecji))

def calcBA(atmcoordi, atmcoordj, atmcoordk):
    """
        Calculate angle ijk of triangle, 
        i.e. the angle between two bonds.
    """
    vecji = atmcoordi - atmcoordj
    vecjk = atmcoordk - atmcoordj
    normji = np.sqrt( np.dot(vecji, vecji) )
    normjk = np.sqrt( np.dot(vecjk, vecjk) )
    dotijk =  np.dot(vecji, vecjk)
    angleijk = 180 - np.arccos( dotijk / (normji * normjk) ) * rad2grad 
    return angleijk

def calcTD(atmcoordi, atmcoordj, atmcoordk, atmcoordl, 
           positive = True):
    """
        Calculate the torsional dihedral between
        the two planes spanned by ijk and jkl. 
        
    """
    atmcoord
    vecij = atmcoordj - atmcoordi
    vecjk = atmcoordk - atmcoordj
    veclk = atmcoordl - atmcoordk
    nrmlijk = np.cross(vecij, vecjk)
    nrmljkl = np.cross(vecjk,veckl)
    ccijkl = np.dot( vecjk, np.cross(nrmlijk, nrmljkl) )
    ddijkl = np.linalg.norm(vecjk) * np.dot(nrmlijk, nrmljkl)
    torsi = np.arctan2(ccijkl, ddijkl) * rad2grad
    if (torsi < 0) and positive:
        torsi += 360
    return torsi

def calcIntrnl(atomcoords, internalType):
    """
        This function simply calls the other
        calc function based on the internalType
        it recieves. 
    """
    intrnl = 0
    if internalType == "bl":
        intrnl = calcBL(atomcoords[0], atomcoords[1])
    elif internalType == "ba":
        intrnl = calcBA(atomcoords[0], atomcoords[1], 
                        atomcoords[2])
    elif internalType == "td":
        intrnl = calcTD(atomcoords[0], atomcoords[1], 
                        atomcoords[2], atomcoords[3])
    return intrnl
