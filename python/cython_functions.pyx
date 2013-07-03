import numpy as np
cimport numpy as np
cimport cython
import scipy.stats as scipy_stats

# type definitions for floats
FDTYPE = np.float64
ctypedef np.float64_t FDTYPE_t

# unsigned short
SDTYPE = np.uint16
ctypedef np.uint16_t SDTYPE_t

# unsigned byte
BDTYPE = np.uint8
ctypedef np.uint8_t BDTYPE_t

@cython.boundscheck(False)

#
# This function takes an input data structure which
# contains label values and a list of unwanted/undesired
# labels. The output data structure contains the input
# labels where the unwanted/undesired labels have been
# replaced by the closest valid label as follows:
#
# For an unwanted label:
#
# 1) all 26 neighbors of a voxel are considered to have
#    the same distance. If 1 or more of these 26 
#    neighbors has a valid label (not 0, nor one of the
#    values in the unwanted list), it will be used
#    to re-label the unwanted label (using the 
#    most occuring label)
#
# Input:  inlabels,       the current set of labels
#         outlables,      a copy of the current set of labels, to be edited
#         unwantedlabels, an array containing the labels numbers that need to be replaced
#
# Output: the function will return the number of labels that could not be replaced,
#         and outlabels will be updated as described above
# 
def substitute_labels_using_26_nearest_neighbors(np.ndarray[SDTYPE_t, ndim=3, mode="c"] inlabels,
                                                 np.ndarray[SDTYPE_t, ndim=3, mode="c"] outlabels,
                                                 np.ndarray[SDTYPE_t, ndim=1] unwantedlabels):
  
  # make sure that the input files have the appropriate types:
  assert inlabels.dtype == SDTYPE and outlabels.dtype == SDTYPE and unwantedlabels.dtype == SDTYPE
  
  # get dimension information from the minc volume
  cdef int nv0 = inlabels.shape[0]
  cdef int nv1 = inlabels.shape[1]
  cdef int nv2 = inlabels.shape[2]
  
  # the number of unwanted labels
  cdef int nlabels = unwantedlabels.shape[0]
  
  # number of unwanted labels that could not be replaced using its 26 neighbors
  cdef int nunresolved = 0

  cdef int v0, v1, v2
  cdef int i0, i1, i2
  
  # from each dimension exclude the first and last element, because these have
  # a number of undefined neighbors: the ones that fall outside of the file
  for v0 in range(1, nv0 - 1):
    for v1 in range(1, nv1 - 1):
      for v2 in range(1, nv2 - 1):
        # first check whether this label value should be replaced
        if( np.rint(inlabels[v0,v1,v2]) in unwantedlabels):
          # create array of possible substitues
          possible_labels = []
          # we can safely look at all 27 voxels in the 3*3 block, because
          # the voxel in question will not be added to the list of 
          # possible labels
          for i0 in (v0-1, v0, v0+1):
            for i1 in (v1-1, v1, v1+1):
              for i2 in (v2-1, v2, v2+1):
                if( not(np.rint(inlabels[i0,i1,i2]) == 0) and not(np.rint(inlabels[i0,i1,i2]) in unwantedlabels)):
                  possible_labels.append(np.rint(inlabels[i0,i1,i2]))
          # potentially we did not find any valid substitues
          if( not(possible_labels) ):
            nunresolved += 1
          else: # we have candidates, simply take the most occuring label number
            outlabels[v0,v1,v2] = np.rint(scipy_stats.mode(possible_labels)[0])
  
  return nunresolved





