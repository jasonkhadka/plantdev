#setup.py file for Quadedge:
from setuptools import setup, Extension
#importing numpy
#import numpy
#getting numpy include directory
#numpy_include = numpy.get_include()

setup(name='quadedge',
    version='0.1',
    ext_modules=[Extension('_quadedge', sources=['edge.cc','list.cc','cell.cc','face.cc','vertex.cc','jacobian.cc','derivative.cc', 'quadedge.i'],
    				 include_dirs=['/usr/local/include','/Users/jasonkhadka/Documents/git/plantdev/nlopt-lib/include'],
    				 library_dirs = ['/usr/local/lib/','/Users/jasonkhadka/Documents/git/plantdev/nlopt-lib/lib/'],
    				 libraries = ["gsl", "nlopt"],
                    swig_opts=['-c++'],
                    )],
    headers=['edge.hh','list.hh','cell.hh','face.hh','vertex.hh','jacobian.hh','derivative.hh']
)
