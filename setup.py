from setuptools import setup
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='STSS',
    version='1.0.1b',
    description='Self-Targeting Spacer Searcher (STSS) identifies self-targeting spacers in CRISPR systems',
    long_description=long_description,
    url='https://github.com/kew222/Self-Targeting-Spacer-Search-tool',
    author='Kyle E. Watters',
    author_email='watters@berkeley.edu',
    license='GNU GPLv3',

    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        ],

    keywords='CRISPR self-targeting',
    py_modules=['STSS','CRISPR_definitions','user_email'],
    install_requires=['requests','biopython'],
    
    )      