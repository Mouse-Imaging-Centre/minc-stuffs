#!/usr/bin/env python3

from setuptools import setup


setup(name="python-stuffs",
      version='0.1.25',
      install_requires=['numpy', 'scipy', 'pyminc'],
      scripts=[f"python/{s}" for s in
              ("TFCE",
               "smooth_vector",
               "measure_xcorr",
               "pmincaverage",
               "minc_label_ops",
               "compute_values_across_segmentation",
               "volumes_from_labels_only",
               "voxel_vote",
               "replace_label_with_nearest_valid_label",
               "rotational_minctracc.py",
               "compute_determinant.py",
               "vtk_meshconvert.py")],
      )
