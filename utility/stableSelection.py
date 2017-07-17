__author__ = 'Haohan Wang'

import numpy as np

def stableSelection(model, X, y, lmbd, lr=None, runs=100, threshold=0.75):
    [n,p] = X.shape
    IND = np.arange(n)
    betaAll = np.zeros(p)
    for i in range(runs):
        ind = np.random.choice(IND, n/2, replace=False)
        X1 = X[ind,:]
        y1 = y[ind]
        model.setLambda(lmbd)
        if lr is not None:
            model.setLearningRate(lr) # learning rate must be set again every time we run it.
        model.fit(X1, y1)
        beta = model.getBeta()
        beta[beta!=0] = 1
        betaAll += beta
    betaAll[betaAll<threshold*n] = 0
    return betaAll