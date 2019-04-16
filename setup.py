
from distutils.core import setup

setup(
      name='cslmm',
      version='0.99',
      author = "Haohan Wang",
      author_email='haohanw@cs.cmu.edu',
      url = "https://github.com/HaohanWang/CS-LMM",
      description = "Discovering Genetic Variants with Weak Associations Guided by Known Variants",
      packages=['models', 'utility'],
      scripts=['cslmm.py'],
    )
