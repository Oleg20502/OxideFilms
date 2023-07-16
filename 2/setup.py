from setuptools import setup, Extension
from Cython.Build import cythonize
import os


dir_path = os.path.dirname(os.path.abspath(__file__))

extension = [Extension('RPDM',
                      sources=[dir_path+'/cfilm.pyx'],
                      include_dirs=['C:/Users/OlegKashurin/anaconda3/include']
            )]

setup(
    name='R-PDM oxide film model class lib',
    ext_modules=cythonize(extension)
)