from distutils.core import setup
import os
from glob import glob
from pybind11.setup_helpers import Pybind11Extension

print("path: ", os.path.abspath(__file__))

# cccrt = Extension(
#     'cccrt',
#     sources=['CCCpy.cpp'],
#     libraries=['boost_python37-mt','boost_numpy37-mt'],
#     extra_compile_args=['-std=c++17'],
#     include_dirs=['..','../hopscotch-map/include'],
# )

cccrt = Pybind11Extension  (
    'cccrt',
    sources=['CCCpy.cpp'],
    depends=['CCC.hpp'],
    extra_compile_args=['-std=c++17', '-O3', '-DPYTHON_BUILD'],
    include_dirs=['..','../eigen-3.4.0', '../EigenRand'],
)

setup(
    name='cccrt',
    version='0.9.1',
    ext_modules=[cccrt])



# setup(
#     name='lutnet',
#     version='1.95',
#     ext_modules=[lutnet])
