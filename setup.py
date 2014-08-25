from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

ext_modules = [Extension("C_Funct", ["C_Funct.pyx"],include_dirs=[numpy.get_include()])]

setup(
  name = 'LSPs',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)




#sudo python setup.py build_ext --inplace