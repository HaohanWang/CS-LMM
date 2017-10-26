
from distutils.core import setup

setup(
      name='cslmm',
      version='0.99',
      author = "Haohan Wang",
      author_email='haohanw@cs.cmu.edu',
      url = "https://github.com/HaohanWang/CS-LMM",
      description = "Discovering Weaker Genetic Associations Guided by Known Associations, with Application to Alcoholism and Alzheimer's Disease Studies",
      packages=['models', 'utility'],
      scripts=['cslmm.py'],
    )
