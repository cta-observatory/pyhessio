===========
pyhessioxxx
===========

CTA Python wrapper for hessio event format that is used in output of simtel_array.

* Code: https://github.com/cta-observatory/pyhessio 

This is a temporaly solution for testing ctapipe with real CTA MC data.
It should not be used once CTA data format DL0 will be accpeted and implemented.

* .. image:: http://img.shields.io/travis/cta-observatory/pyhessio.svg?branch=master
    :target: https://travis-ci.org/cta-observatory/pyhessio
    :alt: Test Status
 
===========
Quick Start
===========

Users of the Anaconda python distribution should follow the instructions for Anaconda python distribution.

Dependencies
------------
    :: 
    python 3 
    setuptools
    numpy

build and install
----------------- 
    ::
    git clone https://github.com/cta-observatory/pyhessio
    cd pyhessio
    python setup.py install 

Installation thanks to anaconda  
________________________________ 
    ::
    git clone https://github.com/cta-observatory/pyhessio
    conda build pyhessio
    conda install --use-local pyhessio
    cd pyhessio

Datasets
____________________________________

To keep the ``pyhessio`` code repo small, we have decided to put the
example files used for unit tests in a separate
repo: https://github.com/cta-observatory/pyhessio-extra ::

    git submodule init
    git submodule update

