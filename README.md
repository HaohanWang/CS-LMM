# CS-LMM (Constrained Sparse multi-locus Linear Mixed Model)

Implementation of CS-LMM in this paper:

    ''Wang H, Aragam B, Lee S, Xing EP, and Wu W. Discovering Weaker Genetic Associations Guided by Known Associations, with Application to Alcoholism and Alzheimer's Disease Studies''

## Introduction

CS-LMM is used to detect the weaker genetic association conditioned on the stronger validated associations.

## File Structure:

* models/ main method for CS-LMM
* utility/ other helper files
* cslmm.py main entry point of using CS-LMM to work with your own data

## An Example Command:

```
python cslmm.py -n data/mice.plink
```
#### Data Support
* CS-LMM currently supports CSV and binary PLINK files.
* Extensions to other data format can be easily implemented through `FileReader` in `utility/dataLoadear`. Feel free to contact us for the support of other data format.

## Installation
You will need to have numpy, scipy and pysnptool installed on your current system.
You can install CS-LMM using pip by doing the following

```
   pip install git+https://github.com/HaohanWang/CS-LMM
```

You can also clone the repository and do a manual install.
```
   git clone https://github.com/HaohanWang/CS-LMM
   python setup.py install
```
## Software with GUI
Software with GUI will be avaliable through [GenAMap](http://genamap.org/)

## Contact
[Haohan Wang](http://www.cs.cmu.edu/~haohanw/)
