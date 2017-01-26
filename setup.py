#setup.py file:
from setuptools import setup, Extension
setup(name='quadedge',
    version='0.1',
    ext_modules=[Extension('_quadedge', sources=['edge.cc','list.cc','cell.cc','face.cc','vertex.cc','jacobian.cc','derivative.cc', 'quadedge.i'],
                    swig_opts=['-c++'],
                    )],
    headers=['edge.hh','list.hh','cell.hh','face.hh','vertex.hh','jacobian.hh','derivative.hh']
)
