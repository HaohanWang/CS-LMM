__author__ = 'Haohan Wang'

import numpy as np
from scipy import stats

class WaldTest:
    def __init__(self, alpha=0.05, fdr=False):
        self.alpha = alpha
        self.fdr = fdr

        
    def tstat(self, beta, var, sigma, q, N, log=False):
        ts = beta / np.sqrt(var * sigma)
        if log:
            ps = 2.0 + (stats.t.logsf(np.abs(ts), N - q))
        else:
            ps = 2.0 * (stats.t.sf(np.abs(ts), N - q))
        return ts, ps

    
    def fit(self, X, y):
        X0 = np.ones(len(y)).reshape(len(y), 1)
        [m, n] = X.shape
        p = []
        for i in range(n):
            Xi = np.hstack([X0 ,X[:, i].reshape(m, 1)])
            XX = np.dot(Xi.T, Xi)
            XX_i = np.linalg.pinv(XX)
            beta = np.dot(np.dot(XX_i, Xi.T), y)
            Uyr = y - np.dot(Xi, beta)
            Q = np.dot(Uyr.T, Uyr)
            sigma = Q * 1.0 / m
            ts, ps = self.tstat(beta[1], XX_i[1, 1], sigma, 1, m)
            if -1e30 < ts < 1e30:
                p.append(ps)
            else:
                p.append(1)
        p = np.array(p)
        self.beta = -np.log(p)

        
    def fdrControl(self):
        tmp = np.exp(-self.beta)
        tmp = sorted(tmp)
        threshold = 1e-8
        n = len(tmp)
        for i in range(n):
            if tmp[i] < (i+1)*self.alpha/n:
                threshold = tmp[i]
        self.beta[self.beta<-np.log(threshold)] = 0

        
    def getBeta(self):
        if not self.fdr:
            self.beta[self.beta < -np.log(self.alpha)] = 0
            return self.beta
        else:
            self.fdrControl()
            return self.beta
