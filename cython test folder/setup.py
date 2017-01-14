from distutils.core import setup
from Cython.Build import cythonize

setup(name="Cy_vicsek", ext_modules=cythonize('Cy_vicsek.pyx'),)
