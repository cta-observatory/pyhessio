from setuptools import setup, Extension
import sys

pyhessio_module = Extension(
    'pyhessio.pyhessioc',
    sources=['pyhessio/src/pyhessio.c',
              'hessioxxx/src/atmprof.c',
              'hessioxxx/src/current.c',
              'hessioxxx/src/dhsort.c',
              'hessioxxx/src/eventio.c',
              'hessioxxx/src/eventio_registry.c',
              'hessioxxx/src/fileopen.c',
              'hessioxxx/src/histogram.c',
              'hessioxxx/src/hconfig.c',
              'hessioxxx/src/moments.c',
              'hessioxxx/src/io_histogram.c',
              'hessioxxx/src/io_history.c',
              'hessioxxx/src/io_simtel.c',
              'hessioxxx/src/io_trgmask.c',
              'hessioxxx/src/straux.c',
              'hessioxxx/src/warning.c',
              'hessioxxx/src/io_hess.c' ],
    include_dirs = ['hessioxxx/include',  '.'],
    define_macros=[('CTA', None), ('CTA_MAX_SC', None)]
  )


NAME = 'pyhessio'
VERSION = '0.5.4'
AUTHOR = 'cta developers'
AUTHOR_EMAIL = 'jean.jacquemier@lapp.in2p3.fr'
URL = 'https://github.com/cta-observatory/pyhessio'
DESCRIPTION = 'CTA Python wrapper for hessio event format that is used in output of simtel_array ..'
LICENSE = 'BSD3'

setup(
    name=NAME,
    packages=['pyhessio'],
    version=VERSION,
    description=DESCRIPTION,
    install_requires=['numpy'],
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    license=LICENSE,
    url=URL,
    classifiers=[
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: BSD License',
    'Operating System :: OS Independent',
    'Programming Language :: C',
    'Programming Language :: Cython',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: Implementation :: CPython',
    'Topic :: Scientific/Engineering :: Astronomy',
    'Development Status :: 3 - Alpha',
    ],

    ext_modules=[pyhessio_module],
)
