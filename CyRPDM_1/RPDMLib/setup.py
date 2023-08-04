from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy
import os


dir_path = os.path.dirname(os.path.abspath(__file__))

extensions = [Extension('RPDM',
                      sources=[dir_path+'/cfilm.pyx'], language='c++',
                      include_dirs=['C:/Users/OlegKashurin/anaconda3/include', numpy.get_include()],
                      extra_compile_args=["-std=c++11", "-O3"]
            )]

setup(
    name='R-PDM oxide film model class lib',
    ext_modules=cythonize(extensions)
)