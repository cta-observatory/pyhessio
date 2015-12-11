===========
pyhessioxxx
===========

CTA Python wrapper for hessio event format that is used in output of simtel_array.

* Code: https://github.com/cta-observatory/pyhessio 

This is a temporaly solution for testing ctapipe with real CTA MC data.
It should not be used once CTA data format DL0 will be accpeted and implemented.

This package is composed of two libraries:

* libhessio: Part of sim_telarray program : https://www.mpi-hd.mpg.de/hfm/~bernlohr/sim_telarray
* pyhessio : libhessio Python wrapper

* .. image:: http://img.shields.io/travis/cta-observatory/pyhessio.svg?branch=add_travis
    :target: https://travis-ci.org/cta-observatory/pyhessio
    :alt: Test Status
 


===========
Quick Start
===========

Anaconda
--------
First install Anaconda python distribution for Python3.4
http://continuum.io/downloads#py34

If you already use anaconda, but with python 2.7, you can create a
second anaconda virtual environment for python3 following the instructions here:
http://continuum.io/blog/anaconda-python-3)::
  
    # only if you already run anaconda: install a new virtualenv for
    # cta development:
    conda create -n cta python=3.4 anaconda
    source activate cta  # to switch to the cta python virtualenv

later you can switch back to your the root environment (or another) by running::
    
    source activate root  
    
Anaconda's distribution doesn't interfere with your local python
distribution if you have one, and can be installed without root
privileges. It contains all the required packages. To "activate"
Anaconda's python, just put it's bin directory in your path: (e.g.
`export PATH=$HOME/anaconda/bin:$PATH`).

After installing anaconda and setting your PATH, run the following to update the packages (for now we have no version restrictions, so the latest ones usually work)::

    conda update --all

CMake
-----
You need to install CMake to build C libhessio library

https://cmake.org


build and install
-----------------
Installation thanks to anaconda:
________________________________
Note: if you do not used a local installation of Anaconda, you need root access on your system::

    git clone https://github.com/cta-observatory/pyhessio
    conda build pyhessio
    conda install --use-local pyhessio
    cd pyhessio
    python setup.py develop

Local installation without anaconda:
____________________________________

Note: Use this procedure if you do not have an local installation of anaconda or no root access on your system::

    git clone https://github.com/cta-observatory/pyhessio
    cd pyhessio
    mkdir build_cmake
    cd build_cmake
    cmake .. -DCMAKE_INSTALL_PREFIX=/home1/jacquem/.local
    make install
    cd ..
    python setup.py install --user --hessio_prefix=/home1/jacquem/.local/lib
    python setup.py develop
