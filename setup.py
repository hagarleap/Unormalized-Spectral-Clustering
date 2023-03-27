from setuptools import Extension, setup, find_packages

setup(name='spkmeans_capi',
     version='1.0',
     description='Python wrapper for Spectral Kmeans C extension',
     packages=find_packages(),
     ext_modules=[Extension('mykmeanssp', sources=['spkmeansmodule.c', 'spkmeans.c'])])

