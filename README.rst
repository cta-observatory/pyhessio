========
pyhessio
========

CTA Python wrapper for hessio event format that is used in output of simtel_array.

* Code: https://github.com/cta-observatory/pyhessio 
* Docs: https://cta-observatory.github.io/pyhessio/

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

Installation for an anaconda env 
________________________________ 

Binary installation:

::

    conda install -c cta-observatory pyhessio

Building from source:

::

    conda create -n cta python=3.5
    source activate cta
    git clone https://github.com/cta-observatory/pyhessio
    cd pyhessio  
    python setup.py install   

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


Using pyhessio
--------------

Warnings: 

1. only ONE simtelarray file could be open at once.

2. This wrapper uses ugly global static variables for development simplicity reason,
and because it is a temporally solution until CTA MC data format exists.

3. PROD3 MC use udge amount of memory, and memory is only free when close_file
is executed. To force close_file execution, we force to use context manager
to instantiate object:

.. code-block:: python

  with pyhessio.open('pyhessio-extra/datasets/gamma_test.simtel.gz') as f:
      f.fill_next_event()

