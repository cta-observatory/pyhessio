from setuptools import setup, Extension

pyhessio_module = Extension(
    'pyhessio.pyhessioc',
    sources=['pyhessio/src/pyhessio.c'],
    include_dirs = ['hessioxxx/include',  '.'],
    libraries=['hessio'],
    define_macros=[('CTA', None), ('CTA_MAX_SC', None)]
)

NAME = 'pyhessio'
VERSION = '0.2'
AUTHOR = 'cta developers'
AUTHOR_EMAIL = 'jean.jacquemier@lapp.in2p3.fr'
URL = 'https://github.com/jacquemier/pyhessio'
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
    ext_modules=[pyhessio_module],
)
