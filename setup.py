#!/usr/bin/env python

from setuptools import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

# setuptools DWIM monkey-patch madness
# http://mail.python.org/pipermail/distutils-sig/2007-September/thread.html#8204
import sys
if 'setuptools.extension' in sys.modules:
    m = sys.modules['setuptools.extension']
    m.Extension.__dict__ = m._Extension.__dict__

setup(name="python-stuffs",
      version='0.1.13',
      scripts=["python/TFCE",
               "python/smooth_vector",
               "python/measure_xcorr",
               "python/pmincaverage",
               "python/minc_label_ops",
               "python/compute_values_across_segmentation",
               "python/volumes_from_labels_only",
               "python/voxel_vote",
               "python/replace_label_with_nearest_valid_label",
               "python/rotational_minctracc.py",
               "minc-scripts-jscholz/vtk_meshconvert.py"],
      cmdclass = {'build_ext': build_ext},
      ext_modules = [Extension("cython_functions", ["python/cython_functions.pyx"], include_dirs=[numpy.get_include()])]
      )
