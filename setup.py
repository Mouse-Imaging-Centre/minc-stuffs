#!/usr/bin/env python

from setuptools import setup


setup(name="python-stuffs",
      version='0.1.22',
      install_requires=['numpy', 'scipy', 'pyminc'],
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
      )
