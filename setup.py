#setup.py file:
from setuptools import setup, Extension
setup(name='quadedge',
    version='0.1',
    ext_modules=[Extension('_quadedge', sources=['array.cc','edge.cc','list.cc','cell.cc','face.cc','obj.cc','vertex.cc','jacobian.cc', 'quadedge.i'],
                    swig_opts=['-c++'],
                    )],
    headers=['array.hh','edge.hh','list.hh','cell.hh','face.hh','obj.hh','vertex.hh','jacobian.hh']
)
