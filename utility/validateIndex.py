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
        print 'There is no validated markers available, first conducting wald test to search most significant SNPs'
        model = WaldTest()
        model.fit(X, y)
        beta = model.getBeta()
        ind = np.where(beta!=0)[0]
        print '\tWald testing identified ', len(ind), 'significantly associated SNPs'
        if len(ind)>5:
            print '\tLRVA will only use the most significant five ones'
            ind = np.argsort(beta)[-5:]
        return ind
