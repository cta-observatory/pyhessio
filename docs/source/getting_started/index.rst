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
