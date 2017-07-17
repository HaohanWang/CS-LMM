__author__ = 'Haohan Wang'
import sys
import numpy as np
from models.WaldTest import WaldTest

def generateValidatedIndex(fileName, Xname, X, y):
    if fileName is not None:
        try:
            text = [line.strip() for line in open(fileName)]
        except:
            print "ERROR: cannot open the validated markers file, at ", fileName
            sys.exit()
        ind = []
        for i in range(Xname):
            if Xname[i] in text:
                ind.append(i)
        return np.array(ind)
    else:
        model = WaldTest()
        model.fit(X, y)
        beta = model.getBeta()
        return np.where(beta!=0)[0]