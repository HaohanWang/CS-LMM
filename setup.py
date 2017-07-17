
from distutils.core import setup

setup(
      name='lrva',
      version='0.99',
      author = "Haohan Wang",
      author_email='haohanw@andrew.cmu.edu',
      url = "https://github.com/HaohanWang/LRVA",
      description = 'Discovering Weaker Genetic Associations with Validated Association, with Studies of Alzheimer\'s Disease and Drug Abuse Disorder',
      packages=['models', 'utility'],
      scripts=['lrva.py'],
    )
