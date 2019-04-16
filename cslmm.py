__author__ = 'Haohan Wang'

# Main file for usage of (CS-LMM) Constrained Sparse multi-locus Linear Mixed Model
# Cite information:
# Wang H, Aragam B, Lee S, Xing EP, and Wu W.
# Discovering Weaker Genetic Associations Guided by Known Associations, with Application to Alcoholism and Alzheimer's Disease Studies
#

def printOutHead(): out.write("\t".join(["RANK", "SNP_ID", "EFFECT_SIZE_ABS"]) + "\n")


def outputResult(rank, id, beta):
    out.write("\t".join([str(x) for x in [rank, id, beta]]) + "\n")


from optparse import OptionParser, OptionGroup

usage = """usage: %prog [options] -n fileName
This program provides the basic usage to CS-LMM, e.g:
python cslmm.py -n data/mice.plink
	    """
parser = OptionParser(usage=usage)

dataGroup = OptionGroup(parser, "Data Options")
modelGroup = OptionGroup(parser, "Model Options")

## data options
dataGroup.add_option("-t", dest='fileType', default='plink', help="choices of input file type")
dataGroup.add_option("-n", dest='fileName', help="name of the input file")
dataGroup.add_option("-v", dest='fileValidated', help="list of the validated markers")

## model options
modelGroup.add_option("--lambda", dest="lmbd", default=None,
                      help="the weight of the penalizer. If neither lambda or snum is given, cross validation will be run.")
modelGroup.add_option("--snum", dest="snum", default=None,
                      help="the number of targeted variables the model selects. If neither lambda or snum is given, cross validation will be run.")
modelGroup.add_option("-s", action='store_true', dest="stable", default=False, help="Stability selection")
modelGroup.add_option('-q', action='store_true', dest='quiet', default=False, help='Run in quiet mode')
modelGroup.add_option('-m', action='store_true', dest='missing', default=False, help='Run without missing genotype imputation')

## advanced options
parser.add_option_group(dataGroup)
parser.add_option_group(modelGroup)

(options, args) = parser.parse_args()


def cross_val_score(clf, X, y, cv=5, learningRate=1):
    scores = []
    [n, p] = X.shape
    b = n/cv
    for i in range(cv):
        ind = np.arange(b) + b*i
        Xtr = np.delete(X, ind, axis=0)
        ytr = np.delete(y, ind, axis=0)
        Xte = X[ind,:]
        yte = y[ind]
        clf.setLearningRate(learningRate)
        clf.fit(Xtr, ytr)
        ypr = clf.predict(Xte)
        s = np.mean(np.square(ypr-yte))
        scores.append(s)
    return scores


def crossValidation(clf, X, y, learningRate):
    minError = np.inf
    minLam = 0
    for i in range(-10, 10):
        lam = np.power(10., i)
        clf.setLambda(lam)
        scores = cross_val_score(clf, X, y, cv=5, learningRate=learningRate)
        score = np.mean(np.abs(scores))
        print lam, score
        if score < minError:
            minError = score
            minLam = lam
    print minLam
    clf.setLambda(minLam)
    clf.setLearningRate(learningRate)
    clf.fit(X, y)
    beta = clf.getBeta()
    return beta


import sys
import numpy as np
from utility.dataLoader import FileReader
from utility.stableSelection import stableSelection
from utility.validateIndex import generateValidatedIndex
from models.CSLMM import CSLMM

fileType = 0
IN = None

if len(args) != 0:
    parser.print_help()
    sys.exit()

outFile = options.fileName + '.output'

print 'Running ... '

reader = FileReader(fileName=options.fileName, fileType=options.fileType, imputation=(not options.missing))
X, Y, Xname = reader.readFiles()

model = CSLMM()

Ind = generateValidatedIndex(options.fileValidated, Xname, X, Y)
model.setKnownInd(Ind)
model.setUp(X, Y)

learningRate = 1e-6
betaM = None

print 'Computation starts ... '

if options.snum is not None: # select for a fix number of variable
    snum = float(options.snum)
    min_lambda_default = 1e-50
    max_lambda_default = 1e50
    patience = 50
    minFactor = 0.75
    maxFactor = 1.25
    min_lambda = min_lambda_default
    max_lambda = max_lambda_default
    iteration = 0
    while min_lambda < max_lambda and iteration < patience:
        iteration += 1
        lmbd = np.exp((np.log(min_lambda) + np.log(max_lambda)) / 2.0)
        if not options.stable:
            if not options.quiet:
                print "\tIter:{}\tlambda:{}".format(iteration, lmbd),
            model.setLambda(lmbd)
            model.setLearningRate(learningRate)  # learning rate must be set again every time we run it.
            model.fit(X, Y)
            beta = model.getBeta()
        else:
            beta = stableSelection(model=model, X=X, y=Y, lmbd=lmbd, lr=learningRate)

        c = len(np.where(np.abs(beta) > 0)[0])  # we choose regularizers based on the number of non-zeros it reports
        if not options.quiet:
            print 'Selected ', c, 'snps',
        if c < minFactor*snum:  # Regularizer too strong
            max_lambda = lmbd
            if not options.quiet:
                print '\tRegularization is too strong, shrink lambda'
        elif c > maxFactor*snum:  # Regularizer too weak
            min_lambda = lmbd
            betaM = beta
            if not options.quiet:
                print '\tRegularization is too weak, enlarge lambda'
        else:
            betaM = beta
            break
elif options.lmbd is not None:
    lmbd = float(options.lmbd)
    if not options.stable:
        model.setLambda(lmbd)
        model.setLearningRate(learningRate)  # learning rate must be set again every time we run it.
        model.fit(X, Y)
        betaM = model.getBeta()
    else:
        betaM = stableSelection(model=model, X=X, y=Y, lmbd=lmbd, lr=learningRate)
else:
    betaM = crossValidation(model, X, Y, learningRate)


ind = np.where(betaM != 0)[0]
bs = betaM[ind].tolist()
xname = []
for i in ind:
    xname.append(i)

beta_name = zip(betaM, Xname)
bn = sorted(beta_name)
bn.reverse()

out = open(outFile, 'w')
printOutHead()

for i in range(len(bn)):
    outputResult(i + 1, bn[i][1], bn[i][0])

out.close()

print '\nComputation ends normally, check the output file at ', outFile
