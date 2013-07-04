/****
 * Needs two inputs: a structure segmentation map and
 * a map with counts. Will then determine the total number
 * of counts in a structure by adding up all counts for
 * each of the labels/segmentations.
 ****/

#include <minc2.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char **argv) {
  mihandle_t    hvol_structures, hvol_counts;
  midimhandle_t dimensions_structures[3], dimensions_counts[3];
  misize_t      size_structures[3], size_counts[3];
  unsigned int  start_structures[3], count_structures[3];
  misize_t      start[3], count[3];
  double        structure_counts[999999], voxel_separations[3];
  int           i;
  double        *counts;
  double        *structures;
  double        volume;

  /* initialization of the start coordinates and dimensions counts */
  start[0] = 0;
  start[1] = 0;
  start[2] = 0;

  count[0] = 10;
  count[1] = 10;
  count[2] = 10;

  if (argc != 3) {
    fprintf(stderr, "Needs two inputs: a structure segmentation map and\na map with counts. Will then determine the total number\nof counts in a structure by adding up all counts for\neach of the labels/segmentations.\n\nThe output is a list of all label numbers and number of counts for\nthis label separated by a comma\n\nUsage: compute_counts_for_labels structures.mnc counts.mnc\n");
    return(1);
  }

  /* read the structure volume */
  if (miopen_volume(argv[1], MI2_OPEN_READ, &hvol_structures) != MI_NOERROR) {
    fprintf(stderr, "Error opening volume: %s\n", argv[1]);
    return(1);
  }

  /* read the counts volume */
  if (miopen_volume(argv[2], MI2_OPEN_READ, &hvol_counts) != MI_NOERROR) {
    fprintf(stderr, "Error opening volume: %s\n", argv[2]);
    return(1);
  }

  /* get dimensions */
  miget_volume_dimensions(hvol_structures, MI_DIMCLASS_SPATIAL, MI_DIMATTR_ALL,
                          MI_DIMORDER_FILE, 3, dimensions_structures);
  miget_volume_dimensions(hvol_counts, MI_DIMCLASS_SPATIAL, MI_DIMATTR_ALL,
                          MI_DIMORDER_FILE, 3, dimensions_counts);
  miget_dimension_sizes(dimensions_structures, 3, size_structures);
  miget_dimension_sizes(dimensions_counts, 3, size_counts);

  miget_dimension_separations(dimensions_structures, MI_ORDER_FILE, 
                              3, voxel_separations);

  /* ensure that both volumes have the same dimension lenghts */
  for (i=0; i < 3; i++) {
    if (size_structures[i] != size_counts[i]) {
      fprintf(stderr, "Dimensions of two volumes must be the same.\n");
    }
  }

  count[0] = size_structures[0];
  count[1] = size_structures[1];
  count[2] = size_structures[2];

  /* Allocating memory */
  counts = (double *) malloc(size_counts[0] *
                             size_counts[1] *
                             size_counts[2] * sizeof(double));

  structures = (double *) malloc(size_structures[0] * 
                                 size_structures[1] *
                                 size_structures[2] * sizeof(double));

  /* Getting hyperslab of the counts file (load full file into memory) */
  if (miget_real_value_hyperslab(hvol_counts, 
                                 MI_TYPE_DOUBLE,
                                 start, 
                                 count,
                                 counts) < 0) {
    fprintf(stderr, "Could not get hyperslab\n");
  }

  /* Getting hyperslab of the segmentation/labels file for the structures (load full file into memory) */
  if (miget_real_value_hyperslab(hvol_structures, 
                                 MI_TYPE_DOUBLE,
                                 start, 
                                 count,
                                 structures) < 0) {
    fprintf(stderr, "Could not get hyperslab\n");
  }

  /* set the count for all structures to 0 */
  for (i=0; i < 999998; i++) {
    structure_counts[i] = 0;
  }

  /* the structure labels/segmentations are assumed to be integer values,
     however, often they are not stored as pure integers. Instead the label
     2 can be stored as 1.98 or 2.02 for instance. This is why we add 0.5
     to the label value and cast it as integer. This way, both 1.98 and 2.02
     will be cast to 2. */
  for (i=0; i < count[0] * count[1] * count[2]; i++) {
    structure_counts[(int)(structures[i] + 0.5)] += counts[i];
  }

  /* labels often are not sequential, e.g., defined label numbers
     could be: 2,4,8,13. To make sure we don't output information 
     about non-existing labels, the counts are only printed for 
     labels that have 1 count in them */
  for (i=0; i < 999998; i++) {
    if (structure_counts[i] != 0) {
      printf("%i, %f\n", i, structure_counts[i]);
    }
  }

  return(0);
}



