from distutils.core import setup
import os
from glob import glob
from pybind11.setup_helpers import Pybind11Extension

print("src path: ", os.path.abspath(__file__))


cccrt = Pybind11Extension  (
    'cccrt',
    sources=['CCCpy.cpp'],
    depends=['CCC.hpp'],
    extra_compile_args=['-std=c++17', '-O3', '-DPYTHON_BUILD'],
    include_dirs=['.','./eigen-3.4.0', './EigenRand'],
)

setup(
    name='cccrt',
    version='0.9.3',
    ext_modules=[cccrt],
    )


