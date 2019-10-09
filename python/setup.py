from distutils.core import setup, Extension
import os

print("path: ", os.path.abspath(__file__))

cccrt = Extension(
    'cccrt',
    sources=['CCCpy.cpp'],
    libraries=['boost_python37-mt','boost_numpy37-mt'],
    extra_compile_args=['-std=c++17'],
    include_dirs=['..'],
)

setup(
    name='cccrt',
    version='0.22',
    ext_modules=[cccrt])
