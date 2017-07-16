__author__ = 'Haohan Wang'

import numpy as np
from numpy import linalg


class Lasso:
    def __init__(self, lam=1., lr=1., tol=1e-5, logistic=False, weighted=False):
        self.lam = lam
        self.lr = lr
        self.tol = tol
        self.decay = 0.5
        self.maxIter = 500
        self.logistic = logistic
        self.weighted = weighted

    def setLambda(self, lam):
        self.lam = lam

    def setLogisticFlag(self, logistic):
        self.logistic = logistic

    def setWeightedFlag(self, weighted):
        self.weighted = weighted

    def setLearningRate(self, lr):
        self.lr = lr

    def setMaxIter(self, a):
        self.maxIter = a

    def setTol(self, t):
        self.tol = t

    def fit(self, X, y):
        if self.logistic:
            if self.weighted:
                self.weights = np.ones_like(y)
                c1 = len(np.where(y == 1)[0])
                c0 = len(np.where(y == 0)[0])
                w1 = float(c1) / (c1 + c0)
                w0 = 1 - w1
                self.weights[y == 1] = w0
                self.weights[y == 0] = w1

        X0 = np.ones(len(y)).reshape(len(y), 1)
        X = np.hstack([X, X0])
        shp = X.shape
        self.beta = np.zeros([shp[1], 1])
        resi_prev = np.inf
        resi = self.cost(X, y)
        step = 0
        while np.abs(resi_prev - resi) > self.tol and step < self.maxIter:
            keepRunning = True
            resi_prev = resi
            while keepRunning:
                prev_beta = self.beta
                pg = self.proximal_gradient(X, y)
                self.beta = self.proximal_proj(self.beta - pg * self.lr)
                keepRunning = self.stopCheck(prev_beta, self.beta, pg, X, y)
                if keepRunning:
                    self.lr = self.decay * self.lr
            step += 1
            resi = self.cost(X, y)
        return self.beta

    def cost(self, X, y):
        if self.logistic:
            if not self.weighted:
                return 0.5 * np.sum(
                    np.square(y - 1. / (1 + np.exp(-np.dot(X, self.beta)))).transpose()) + self.lam * linalg.norm(
                    self.beta, ord=1)
            else:
                return 0.5 * np.sum(np.square(
                    y - 1. / (1 + np.exp(-np.dot(X, self.beta)))).transpose() * self.weights) + self.lam * linalg.norm(
                    self.beta, ord=1)
        else:
            return 0.5 * np.sum(np.square(y - np.dot(X, self.beta)).transpose()) + self.lam * linalg.norm(
                self.beta, ord=1)

    def proximal_gradient(self, X, y):
        if self.logistic:
            if not self.weighted:
                return -np.dot(X.transpose(), (y.reshape((y.shape[0], 1)) - 1. / (1 + np.exp(-np.dot(X, self.beta)))))
            else:
                return -np.dot(X.transpose(), (
                    y.reshape((y.shape[0], 1)) - 1. / (1 + np.exp(-np.dot(X, self.beta)))) * self.weights.reshape(
                    self.weights.shape[0], 1))
        else:
            return -np.dot(X.transpose(), (y.reshape((y.shape[0], 1)) - (np.dot(X, self.beta))))

    def proximal_proj(self, B):
        t = self.lam * self.lr
        zer = np.zeros_like(B)
        result = np.maximum(zer, B - t) - np.maximum(zer, -B - t)
        return result

    def predict(self, X):
        X0 = np.ones(X.shape[0]).reshape(X.shape[0], 1)
        X = np.hstack([X, X0])
        return np.dot(X, self.beta)

    def getBeta(self):
        self.beta = self.beta.reshape(self.beta.shape[0])
        return self.beta[:-1]

    def stopCheck(self, prev, new, pg, X, y):
        if np.square(linalg.norm((y - (np.dot(X, new))))) <= \
                                np.square(linalg.norm((y - (np.dot(X, prev))))) + np.dot(pg.transpose(), (
                                    new - prev)) + 0.5 * self.lam * np.square(linalg.norm(prev - new)):
            return False
        else:
            return True
