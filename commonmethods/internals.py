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

class molGraph(object):
    # Simple molecular graph class, where vertices are
    # atoms and edges are bonds 
    def __init__(self,graphinput=None):
        if graphinput == None:
            graphinput = {}
        self.__graph_dict = graphinput
    
    def add_vertex(self,vertex):
        if vertex not in self.__graph_dict:
            self.__graph_dict[vertex] = []  
            return True
        else:
            return False
    
    def add_edge(self,edge):
        (vertex1, vertex2) = (edge[0], edge[1])
        if vertex1 in self.__graph_dict:
            self.__graph_dict[vertex1].append(vertex2)
        else:
            self.__graph_dict[vertex1] = [vertex2]

    def get_vertices(self):
        return list(self.__graph_dict.keys())

    def get_edges(self):
        edges = []
        for node in self.__graph_dict:
            for neighbour in self.__graph_dict[node]:
                edges.append({node, neighbour})

        return edges

    def get_neighbours(self,vertex):
        return self.__graph_dict[vertex]

    def get_dict(self):
        return self.__graph_dict

    def depth_first_search(self,vertex,paths,discovered,curPathLength,
                           maxPathLength=4):
        """ This is a modified depth first search method that determines 
            all simple paths of a certain length. In the usual DFS one
            discovers vertices only once, while here multiple discoveries
            are allowed. The list paths is passed by reference and contains 
            the output of the function.
            Inputs: vertex (str) and maximum pathlength (int)
            Outputs: list of paths with maximum pathlengt (list of strs)    """
        curPathLength += 1   # each recursion, increase current path length 
        neighbours = self.get_neighbours(vertex)
        # only recurse if maximum pathlength is not yet reached
        if curPathLength < maxPathLength: 
            for neighbour in neighbours:
                # only recurse if the vertex has not yet been discovered
                if neighbour not in discovered:
                    # add current neighbour to the discovered list
                    discovered.append(neighbour)
                    # recurse 
                    test = self.depth_first_search(neighbour,paths,discovered,
                                            curPathLength, maxPathLength)
                # only if the bottom is actually reached will test a value
                    if test != None:
                        tmp = test[:]
                        paths.append(tmp)
                    discovered.remove(neighbour)
        elif curPathLength == maxPathLength:
            # bottom is reached 
            return discovered

class redundantInternals(molGraph):
    """ This class is inteded for conversion of Cartesian coordinates 
        to redundant internals. Basically the cartesians are read from
        some file (right now only Turbomole is supported -> change in
        future) and the relevant bond lengths are calculated first. 
        From this the molecular graph molGraph class is constucted.
        Based on this connectivity all bond angles and torsional 
        dihedrals are determined. """
    def __init__(self, numAtoms, fname=None, fileType=None, inpMolGraph=None,
                 mollist=None):
        if fname != None and mollist == None and inpMolGraph==None:
            self.setupGraph = True
            self.numAtoms = numAtoms
            assert not(fileType == None)
            self.fileType = fileType 
            self.atomNames, self.atomCoords = self.getCartesians(fname)
            self.distMat = self.initDistanceMatrix()
            self.BLConnectivity, self.bondLengths = self.initBondLengths()
        elif fname != None and mollist == None and inpMolGraph!=None:
            self.setupGraph = False 
            self.numAtoms = numAtoms
            self.atomNames, self.atomCoords = self.getCartesians(fname)
            self.distMat = self.initDistanceMatrix()
        elif fname == None and mollist != None:
            self.atomNames, self.atomCoords = mollist 
            self.numAtoms = numAtoms
            self.distMat = self.initDistanceMatrix()

        if inpMolGraph == None:
            molGraph.__init__(self)
            for indBL, BL in enumerate(self.BLConnectivity):
                self.add_vertex(indBL)
                self.add_edge([indBL , BL])
                self.add_edge([BL , indBL])
        else:
            molGraph.__init__(self,graphinput=inpMolGraph)
            if self.atomCoords != None:
                self.BLConnectivity, self.bondLengths = self.initBondLengths()

    def getCartesians(self, fileName):
        with open(fileName, "r") as fileLines:
            atomName = []
            atomCoord = []
            if self.fileType == "xyz":  
                nrLines = -1 
            elif self.fileType == "dat":
                nrLines = -2
            elif self.fileType == "tm":
                nrLines = -1
            for fileLine in fileLines:
                nrLines += 1
                if (self.fileType == "xyz") or (self.fileType == "dat"): 
                    if nrLines >= 1:
                        atomName.append(fileLine.strip().split()[0] + str(nrLines))
                        atomCoord.append([float(fileLine.strip().split()[i])*b2A
                                          for i in range(1,4)])
                elif (self.fileType == "tm"):
                    atomName.append(fileLine.strip().split()[3] + str(nrLines))
                    atomCoord.append([float(fileLine.strip().split()[i])*b2A
                                      for i in range(3)])
                    
                if nrLines == self.numAtoms:
                    break
        return atomName, atomCoord 
    
    def initDistanceMatrix(self):
        distMat = np.zeros((self.numAtoms,self.numAtoms))
        for i in range(self.numAtoms):
            for j in range(i+1,self.numAtoms):
                posi = np.array(self.atomCoords[i])
                posj = np.array(self.atomCoords[j])
                vecji = posi - posj 
                distji = np.sqrt(np.dot(vecji,vecji))
                distMat[i,j] = distji
                distMat[j,i] = distji

        return distMat

    def initBondLengths(self):
        print("Bond lengths:")
        bondLengths = []
        distconnect = np.zeros(self.numAtoms).astype(int)
        dist = np.zeros(self.numAtoms)
        if self.setupGraph:
            for i in range(1,self.numAtoms):
                indMin, distMin = (np.argmin(self.distMat[i,:i]),
                                   np.amin(self.distMat[i,:i]))
                atomInit = self.atomNames[indMin]
                atomFinal = self.atomNames[i]
                #print('{strAtomi:3} - {strAtomj:4}: {strDistji:.3f} pm'.format(
                #      strAtomi = atomInit, strAtomj = atomFinal,
                #      strDistji = distMin*100))
                distconnect[i] = indMin 
                dist[i] = distMin
        else:
            for i in range(1,self.numAtoms):
                bondPartner = self.get_neighbours(i)[0]
                bondLength = self.distMat[i, self.get_neighbours(i)[0]]
                distconnect[i] = bondPartner
                dist[i] = bondLength 
                atomInit = self.atomNames[bondPartner]
                atomFinal = self.atomNames[i]
                #print('{strAtomi:3} - {strAtomj:4}: {strDistji:.3f} pm'.format(
                #      strAtomi = atomInit, strAtomj = atomFinal,
                #      strDistji = bondLength*100))
        
        return distconnect, dist

    def getBondPartner(self): 
        for i in np.arange(1, self.BLConnectivity.size):
            print(self.atomNames[self.BLConnectivity[i]], self.atomNames[i])
