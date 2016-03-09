========
pyhessio
========

CTA Python wrapper for hessio event format that is used in output of simtel_array.

* Code: https://github.com/cta-observatory/pyhessio 

This is a temporaly solution for testing ctapipe with real CTA MC data.

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

    python3 
    setuptools
    numpy

build and install
----------------- 

Installation without anaconda  
________________________________ 
::
    git clone https://github.com/cta-observatory/pyhessio
    cd pyhessio
    python setup.py install 

Installation for an anaconda env 
________________________________ 

::

    conda create -n cta python=3.4
    source activate cta
    git clone https://github.com/cta-observatory/pyhessio
    conda build pyhessio
    conda install --use-local pyhessio

Datasets
____________________________________

To keep the ``pyhessio`` code repo small, we have decided to put the
example files used for unit tests in a separate
repo: https://github.com/cta-observatory/pyhessio-extra ::

    git submodule init
    git submodule update

For pyhessio  developers
________________________

::

    python setup.py develop
