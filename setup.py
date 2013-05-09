#! /usr/bin/env python

# System imports
from distutils.core import *
from distutils      import sysconfig
from distutils.extension import Extension

# Third-party modules - we depend on numpy for everything
import numpy

# Obtain the numpy include directory.
numpy_include = numpy.get_include()
_phasefit = Extension("_phasefit",
        ["phasefit_wrap.cxx",
          "phasefit.cc"],
        include_dirs = [numpy_include],
        )
_maskarray = Extension("_maskarray",
        ["maskarray_wrap.cxx",
          "maskarray.cc"],
        include_dirs = [numpy_include],
        )
_transitmask = Extension("_transitmask",
        ["transitmask_wrap.cxx",
          "transitmask.cc"],
        include_dirs = [numpy_include],
        )
_fitmap = Extension("_fitmap",
        ["fitmap_wrap.cxx",
          "fitmap.cc"],
        include_dirs = [numpy_include],
        )

_blssimple = Extension("_blssimple",
        ["blssimple_wrap.cxx",
          "blssimple.cc"],
        include_dirs = [numpy_include],
        )
setup(name= "Outlier",
	description = "phasefold to deal with outliers",
	author      = "Xu Huang",
    author_email = "chelsea@astro.princeton.edu",
    url = "",
    version     = "1.0.0",
	py_modules  = ["phasefit","maskarray","transitmask","fitmap","blssimple"],
    ext_modules = [_phasefit,_maskarray,_transitmask,_fitmap,_blssimple])

