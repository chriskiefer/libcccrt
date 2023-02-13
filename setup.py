from distutils.core import setup
import os
from glob import glob
from pybind11.setup_helpers import Pybind11Extension, build_ext

# print("src path: ", os.path.abspath(__file__))


cccrt = Pybind11Extension  (
    'cccrt',
    sources=['CCCpy.cpp'],
    depends=['CCC.hpp'],
    extra_compile_args=['-O3', '-DPYTHON_BUILD', '$(python3.10-config --includes)', '$(python3.10 -m pybind11 --includes)'],
    include_dirs=['.','./eigen-3.4.0', './EigenRand'],
)

setup(
    name='cccrt',
    version='0.9.3',
    cmdclass={"build_ext": build_ext},
    ext_modules=[cccrt],
    )


