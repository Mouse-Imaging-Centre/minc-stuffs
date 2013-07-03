#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext


setup(name="python-stuffs",
      version='0.1.1',
      scripts=["python/TFCE", 
	       "python/smooth_vector", 
	       "python/measure_xcorr", 
               "python/pmincaverage", 
               "python/minc_label_ops", 
               "python/compute_values_across_segmentation", 
               "python/volumes_from_labels_only", 
               "python/voxel_vote"],
      cmdclass = {'build_ext': build_ext},
      ext_modules = [Extension("cython_functions", ["python/cython_functions.pyx"])]
      )
