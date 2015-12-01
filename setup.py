from setuptools import setup, Extension
import sys

hessio_prefix=None

# Search for hessio_prefix for local install
if any("hessio_prefix" in s for s in sys.argv):
  hessio_prefix = [s for s in sys.argv if 'hessio_prefix' in s][0]
  sys.argv.remove(hessio_prefix)
  hessio_prefix = hessio_prefix.split('=')[1]

  pyhessio_module = Extension(
      'pyhessio.pyhessioc',
      sources=['pyhessio/src/pyhessio.c'],
      include_dirs = ['hessioxxx/include',  '.'],
      libraries=['hessio'],
      library_dirs=[hessio_prefix],
      runtime_library_dirs=[hessio_prefix],
      define_macros=[('CTA', None), ('CTA_MAX_SC', None)]
  )

# else standars installation
else:
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
    ext_modules=[pyhessio_module],
)

