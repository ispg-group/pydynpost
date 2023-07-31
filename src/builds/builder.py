#!/usr/bin/env python
import copy
import sympy as sp
import numpy as np
from sympy.utilities.autowrap import ufuncify
import src.traj as traj
import src.builds.vectorize as vectorize
import src.builds.pseudospec as pseudospec

class Builder():
    def calculateS(self, B1, B2):
        enumTrajPairs = traj.constructPairIterator(B1.trajs, B2.trajs)
        trajOverlap = np.zeros((B1.nTraj, B2.nTraj), dtype=complex)
        densityMat = np.zeros((B1.nTraj, B2.nTraj), dtype=complex)
        bundleOverlap = 0.0 + 0.0j
        for iTrajPair, trajPair in enumTrajPairs:
            traj1 = trajPair[0]
            traj2 = trajPair[1]
            trajOverlap[iTrajPair] = self.calculateSprim(traj1,traj2)
            densityMat[iTrajPair] = np.conjugate(traj1.amp)*traj2.amp
            bundleOverlap += trajOverlap[iTrajPair] \
                             * densityMat[iTrajPair]

        return trajOverlap, densityMat, bundleOverlap

    def calculateH(self, B1, B2, trajOverlap, densityMat):
        enumTrajPairs = traj.constructPairIterator(B1.trajs, B2.trajs)
        trajT = 0.0 + 0.0j
        bundleT = 0.0 + 0.0j 
        bundleV = 0.0 + 0.0j
        bundleC = 0.0 + 0.0j
        bundleSDot = 0.0 + 0.0j
        for iTrajPair, trajPair in enumTrajPairs:
            traj1 = trajPair[0]
            traj2 = trajPair[1]
            trajT = self.calculateTprim(
                traj1,traj2,trajOverlap[iTrajPair]
            )
            bundleT += trajT * densityMat[iTrajPair] 
            if traj1.statID == traj2.stateID: 
                trajV = pseudospec.Potential(traj1, traj2, trajOverlap)
                bundleV += trajV * densityMat[iTrajPair] 
            else:
                trajC = pseudospec.Coupling(traj1, traj2, trajOverlap)
                bundleC += trajC * densityMat[iTrajPair] 

            densityDot = self.calculateDensTerm(traj1, traj2)
            overlapDot = self.calculateOverlapTerm(
                traj1, traj2, trajOverlap
            )
            trajSDot = self.calculateSDotprim(
                densityDot,trajOverlap[iTrajPair],  
                overlapDot,densityMat[iTrajPair]
            )
            bundleSDot += trajSDot
             

    def calculateSprim(self, traj1, traj2):
        a_i = traj1.widths
        a_j = traj2.widths
        Q_i = self.flatten(traj1.positions)
        Q_j = self.flatten(traj2.positions)
        P_i = self.flatten(traj1.momenta)
        P_j = self.flatten(traj2.momenta)
        gamma_i = traj1.phase
        gamma_j = traj2.phase
        oneDOverlaps = self.vectorizedMat.S(
            a_i, a_j, Q_i, Q_j, P_i, P_j
        ) 
        fullDOverlap = np.prod(oneDOverlaps) 
        fullDOverlap *= np.exp(1.j*(gamma_j - gamma_i))
        return fullDOverlap

    def calculateTprim(self, traj1, traj2, S_ij):
        a_i = traj1.widths
        a_j = traj2.widths
        Q_i = self.flatten(traj1.positions)
        Q_j = self.flatten(traj2.positions)
        P_i = self.flatten(traj1.momenta)
        P_j = self.flatten(traj2.momenta)
        gamma_i = traj1.phase
        gamma_j = traj2.phase
        oneDKinetic = self.vectorizedMat.T(
            a_i, a_j, Q_i, Q_j, P_i, P_j
        ) 
        fullDKinetic = np.prod(oneDKinetic) * S_ij 
        return fullDKinetic

    
    def calculateOverlapTerm(self, traj1, traj2, S_ij):
        a_i = traj1.widths
        a_j = traj2.widths
        Q_i = self.flatten(traj1.positions)
        Q_j = self.flatten(traj2.positions)
        P_i = self.flatten(traj1.momenta)
        P_j = self.flatten(traj2.momenta)
        gamma_i = traj1.phase
        gamma_j = traj2.phase
        F_j = self.flatten(traj2.momenta)
        oneDdSdQ_j = self.vectorizedMat.dSdQ(
            a_i, a_j, Q_i, Q_j, P_i, P_j
        )
        oneDdSdP_j = self.vectorizedMat.dSdP(
            a_i, a_j, Q_i, Q_j, P_i, P_j
        )
        dSdQ_j = np.prod(oneDdSdQ_j) * S_ij
        dSdP_j = np.prod(oneDdSdP_j) * S_ij
        SDot = np.dot(P_j/M_j,dSdQ_j) \
             + np.dot(F_j, dSdP_j)    \

    def setVectorized(self, vectorizedMat):
        self.vectorizedMat = vectorizedMat

    def flatten(self, inDict): 
        flattenedVec = []
        for value in inDict.values():
            flattenedVec.extend(value.tolist())

        return np.array(flattenedVec)


def generateMatrixElements(builder):
    alpha_i, alpha_j = sp.symbols(
        'alpha_i alpha_j', real=True, positive=True 
    ) 
    p_i, q_i = sp.symbols('p_i q_i', real=True)
    p_j, q_j = sp.symbols('p_j q_j', real=True)
    q = sp.symbols('q', real=True)
    vectorizedMat = vectorize.VectorizedMat(
        alpha_i, alpha_j, q_i, q_j, p_i, p_j, q
    )
    outBuilder = copy.deepcopy(builder)
    outBuilder.setVectorized(vectorizedMat)
    return outBuilder 
    
