__author__ = 'Haohan Wang'

import numpy as np

from Lasso import Lasso

class LRVA:
    def __init__(self, lam=1., lr1=1e-1, lr2=1e-6, tol=1e-5, maxIter=500, snpFile=True, logistic=False, weighted=False):
        self.lam = lam
        self.lr1 = lr1
        self.lr2 = lr2
        self.tol = tol
        self.maxIter = maxIter
        self.decay = 0.5
        self.snpFile = snpFile
        self.logistic = logistic
        self.weighted = weighted

    def setLambda(self, ldb):
        self.lam = ldb

    def setLogisticFlag(self, logistic):
        self.logistic = logistic

    def setWeightedFlag(self, weighted):
        self.weighted = weighted

    def setLearningRate(self, lr):
        self.lr2 = lr

    def setLearningRate1(self, lr):
        self.lr1 = lr

    def setTolerance(self, tol):
        self.tol = tol

    def setMaxIter(self, m):
        self.maxIter = m

    def setKnownInd(self, ind):  # set the known associations with index, 0, 2, 3 etc.
        self.kI = ind

    def setSnpFlag(self, snpFile):
        self.snpFile = snpFile

    def calculateLinearDependentCorrelatedVariables(self, X):
        [m, n] = X.shape
        result = []
        if not self.snpFile:
            for i in range(n):
                if i not in self.kI:
                    X2 = X[:, i]
                    C11 = np.dot(X2.T, X2) * 1.0 / n
                    C21 = np.dot(self.X1.T, X2) * 1.0 / n
                    ii = 1.0 / C11
                    r = np.abs(np.dot(C21, ii))
                    c = len(np.where(r >= 1)[0])
                    if c > 0:
                        result.append(i)
                    else: # if there is no linear dependent relationship, test for correlation
                        for j in range(self.X1.shape[1]):
                            col = self.X1[:,j]
                            cor = np.corrcoef(col, X2)
                            if np.abs(cor[0][1]) > 0.9:
                                result.append(i)
                            break
        else:
            pass
        return result

    def cross_val_score(self, clf, X, y, cv=5):
        scores = []
        [n, p] = X.shape
        b = n/cv
        for i in range(cv):
            ind = np.arange(b) + b*i
            Xtr = np.delete(X, ind, axis=0)
            ytr = np.delete(y, ind, axis=0)
            Xte = X[ind,:]
            yte = y[ind]
            clf.fit(Xtr, ytr)
            ypr = clf.predict(Xte)
            if np.mean(np.abs(ypr)) == 0:
                s = 1e100
            else:
                s = np.mean(np.square(ypr-yte))
            scores.append(s)
        return scores

    def fitBeta(self, X, y):
        self.phase1model = Lasso(lam=0, logistic=self.logistic, weighted=self.weighted)
        self.phase1model.setLearningRate(self.lr1)
        self.phase1model.fit(X, y)
        beta = self.phase1model.getBeta()
        yr = self.phase1model.predict(X)
        return beta, yr

    def setUp(self, X, y):
        self.y = y
        [n, p] = X.shape
        # setup
        self.kIComplementary = []
        self.X1 = X[:, self.kI]

        self.kINone = self.calculateLinearDependentCorrelatedVariables(X)
        for i in range(p):
            if i not in self.kI:
                if i not in self.kINone:
                    self.kIComplementary.append(i)
        self.X2 = X[:, self.kIComplementary]

        # phase one
        self.b1, yr = self.fitBeta(self.X1, y)
        self.c = np.min(np.abs(self.b1))
        if self.logistic:
            y_tmp = y+1e-5
            self.y2 = -np.log(np.abs(1-(y_tmp))/(y_tmp)) - yr
        else:
            self.y2 = y - yr
        self.bias = yr
        self.nkI = []

    def assemble(self):
        p = len(self.kI) + len(self.kIComplementary) + len(self.kINone)
        self.beta = np.zeros([p])
        self.beta[self.kI] = self.b1
        self.beta[self.kIComplementary] = self.b2

    def fit(self, X=None, y=None):
        self.phase2model = Lasso()
        self.phase2model.setLearningRate(self.lr2)
        self.phase2model.setLambda(self.lam)
        self.phase2model.setMaxIter(self.maxIter)
        self.phase2model.setTol(self.tol)
        self.phase2model.fit(self.X2, self.y2)

        self.b2 = self.phase2model.getBeta()

    def getBeta(self):
        self.assemble()
        return self.beta

    def predict(self, X=None):
        Xtmp1 = X[:, self.kI]
        Xtmp2 = X[:, self.kIComplementary]
        return self.phase2model.predict(Xtmp2) + self.phase1model.predict(Xtmp1)
