#!/usr/bin/env python

from distutils.core import setup

setup(name="python-stuffs",
      version='0.1',
      scripts=["python/TFCE", 
	       "python/smooth_vector", 
	       "python/measure_xcorr", 
               "python/pmincaverage", 
               "python/minc_label_ops", 
               "python/average_values_across_segmentation", 
               "python/volumes_from_labels_only", 
               "python/voxel_vote"]
      )
