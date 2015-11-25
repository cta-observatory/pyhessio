#!/bin/sh

cmake . -DCMAKE_INSTALL_PREFIX=$PREFIX
make -j$(getconf _NPROCESSORS_ONLN)
make install

export CPATH=$PREFIX/include/pyhessio
export LIBRARY_PATH=$PREFIX/lib
python setup.py install
python setup.py develop
