from distutils.core import setup, Extension

module1 = Extension(
    'algorithm',
    sources=['./algorithm/algorithmmodule.c', './algorithm/algorithm.c'],
)

setup(name='biosequence',
      version='1.0',
      description='this is a demo package',
      ext_modules=[module1],
      packages=['biosequence'])
