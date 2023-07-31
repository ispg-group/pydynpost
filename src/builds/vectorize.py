#!/usr/bin/env python
import sympy as sp

class OneDGaussian():
    def __init__(self,alpha, q, qc, pc):
        self.alpha = alpha
        self.q = q
        self.qc = qc
        self.pc = pc
        self.oneDGaussian = self.setupOneDGaussian(alpha, q, qc, pc) 

    def setupOneDGaussian(self, alpha, q, qc, pc):
        root = 2 * alpha / sp.pi
        exponent = - alpha * (q - qc)**2 + sp.I * pc* (q - qc)
        oneDGaussian = sp.real_root(root, 4)
        oneDGaussian *= sp.exp(exponent) 
        return oneDGaussian

    def vectorize(self, alphaIn, qcIn, pcIn):
        partialFunc = self.oneDGaussian.subs(
            [(self.alpha,alphaIn), (self.qc,qcIn), (self.pc,pcIn)]
        )
        return ufuncify(self.q, partialFunc, backend='f2py')
    
    def __call__(self):
        return self.oneDGaussian

class TangentOneDGaussian():
    def __init__(self,alpha, q, qc, pc):
        self.alpha = alpha
        self.q = q
        self.qc = qc
        self.pc = pc
        oneDGaussian = OneDGaussian(alpha, q, qc, pc) 
        self.tangentOneDGaussian = diff(self.oneDGaussian,self.q)  

    def vectorize(self, alphaIn, qcIn, pcIn):
        partialFunc = self.tangentOneDGaussian.subs(
            [(self.alpha,alphaIn), (self.qc,qcIn), (self.pc,pcIn)]
        )
        return ufuncify(self.q, partialFunc, backend='f2py')

class VectorizedMat():
    def __init__(self, alpha_i, alpha_j, q_i, q_j, p_i, p_j, q): 
        self.oneDGaussian_i = OneDGaussian(alpha_i, q, q_i, p_i) 
        self.oneDGaussian_j = OneDGaussian(alpha_j, q, q_j, p_j) 
        self._S = self.deriveS() 
        self.S = self.vectorizeS()  
        self._T = self.deriveT()
        self.T = self.vectorizeT()
        self._dSdQ = self.derivedSdQ() 
        self.dSdQ = self.vectorizedSdQ()
        self._dSdP = self.derivedSdP() 
        self.dSdP = self.vectorizedSdP()
    
    def _vectorize(self, mat):
        return ufuncify(
            (self.oneDGaussian_i.alpha, self.oneDGaussian_j.alpha,
            self.oneDGaussian_i.qc, self.oneDGaussian_j.qc,
            self.oneDGaussian_i.pc, self.oneDGaussian_j.pc),
            mat, backend='f2py' 
        ) 

    def deriveS(self): 
        prod = sp.simplify(
            sp.conjugate(self.oneDGaussian_i()) * self.oneDGaussian_j()
        )
        overlap = sp.simplify(
            sp.integrate(prod, (self.oneDGaussian_j.q, -sp.oo, sp.oo))
        )
        prefactorTerms = overlap.args[:-1]
        prefactor = 1
        for prefactorTerm in prefactorTerms:
            prefactor *= prefactorTerm
        exponential = sp.simplify(overlap.args[-1])
        exponentTerms = exponential.args[0].args   
        numerator = 0
        denominator = 1
        for exponentTerm in exponentTerms:
            if exponentTerm.func == sp.Add:
                numerator = exponentTerm 
            else:
                denominator *= exponentTerm
        tmpNumerator = sp.expand(numerator) 
        ai = self.oneDGaussian_i.alpha
        aj = self.oneDGaussian_j.alpha
        aiMaj = ai * aj
        aiMajTerm = aiMaj * tmpNumerator.coeff(aiMaj) 
        tmpNumerator = sp.simplify(tmpNumerator - aiMajTerm)
        aiMajTerm = sp.factor(aiMajTerm)
        aiTerm = sp.collect(tmpNumerator, ai).coeff(ai,1) * ai 
        ajTerm = sp.collect(tmpNumerator, aj).coeff(aj,1) * aj 
        aiPajTerm = aiTerm + ajTerm 
        tmpNumerator = sp.simplify(tmpNumerator - aiPajTerm)
        aiPajTerm = sp.factor(aiPajTerm)
        finalTerm = sp.factor(tmpNumerator) 
        numerator = aiMajTerm + aiPajTerm + finalTerm 
        exponent = numerator * denominator
        overlap = prefactor * sp.exp(exponent)
        return overlap

    def vectorizeS(self):
        return self._vectorize(self._S)

    def deriveT(self):
        d2Sdq2 = sp.diff(self._S, self.oneDGaussian_j.qc, 2) 
        polynomial = sp.cancel(d2Sdq2 / self._S)
        denominator = sp.factor(polynomial.args[0]) 
        numerator = polynomial.args[1]
        reTerm = numerator.coeff(sp.I,0)
        imTerm = numerator.coeff(sp.I,1)

        ai = self.oneDGaussian_i.alpha
        aj = self.oneDGaussian_j.alpha
        ai2Maj2 = ai**2 * aj**2
        ai2Maj2Term = reTerm.coeff(ai2Maj2, 1) * ai2Maj2  
        tmpReTerm = sp.simplify(reTerm - ai2Maj2Term)
        ai2Maj2Term = sp.factor(ai2Maj2Term)
        aiMaj2 = ai * aj**2
        ai2Maj = ai**2 * aj
        aiMaj2Term = tmpReTerm.coeff(aiMaj2, 1) * aiMaj2
        aiMaj2Term += tmpReTerm.coeff(ai2Maj, 1) * ai2Maj
        tmpReTerm = sp.simplify(tmpReTerm - aiMaj2Term)
        aiMaj2Term = sp.factor(aiMaj2Term)
        aiMajTerm = sp.factor(tmpReTerm) 
        reTerm = ai2Maj2Term + aiMaj2Term + aiMajTerm

        aiMaj2Term = imTerm.coeff(aiMaj2, 1) * aiMaj2
        aiMaj2Term += imTerm.coeff(ai2Maj, 1) * ai2Maj
        imTerm = sp.factor(aiMaj2Term)
    
        numerator = reTerm + sp.I * imTerm
        d2Sdq2 = numerator * denominator 
        return -d2Sdq2

        
    def vectorizeT(self):
        return self._vectorize(self._T)

    def derivedSdQ(self):
        dSdQ = sp.diff(self._S, self.oneDGaussian_j.qc, 1) 
        polynomial = sp.cancel(dSdQ / self._S)
        dSdQ = polynomial 
        return dSdQ
    
    def vectorizedSdQ(self):
        return self._vectorize(self._dSdQ)

    def derivedSdP(self):
        dSdP = sp.diff(self._S, self.oneDGaussian_j.pc, 1)
        polynomial = sp.cancel(dSdP / self._S)
        dSdP = polynomial
        return dSdP

    def vectorizedSdP(self):
        return self._vectorize(self._dSdP)
         
