# LRVA (Linear Regression with Validated Association)

Implementation of LRVA in this paper:

    ''Wang H, Aragam B, Lee S, Xing EP, and Wu W. Discovering Weaker Genetic Associations with Validated Association, with Studies of Alzheimer's Disease and Drug Abuse Disorder''

## Introduction

LRVA is used to detect the weaker genetic association conditioned on the stronger validated associations.

## File Structure:

* models/ main method for LRVA
* utility/ other helper files
* lrva.py main entry point of using LRVA to work with your own data

## An Example Command:

```
python lrva.py -n data/mice.plink
```
#### Data Support
* LRVA currently supports CSV and binary PLINK files.
* Extensions to other data format can be easily implemented through `FileReader` in `utility/dataLoadear`. Feel free to contact us for the support of other data format.

## Installation
You will need to have numpy, scipy and pysnptool installed on your current system.
You can install LRVA using pip by doing the following

```
   pip install git+https://github.com/HaohanWang/LRVA
```
This should make the module pylmm available as well as the two scripts pylmmGWAS.py and pylmmKinship.py.

You can also clone the repository and do a manual install.
```
   git clone https://github.com/HaohanWang/LRVA
   python setup.py install
```

## Contact
[Haohan Wang](http://www.cs.cmu.edu/~haohanw/)
