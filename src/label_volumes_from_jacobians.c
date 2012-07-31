/****
 * Needs two inputs: a structure segmentation map and
 * jacobians. Will then estimate the volume of each structure
 * by multiplying each voxel by the voxel's volume times the 
 * jacobian
 ****/

#include <minc2.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char **argv) {
  mihandle_t    hvol_structures, hvol_jacobians;
  midimhandle_t dimensions_structures[3], dimensions_jacobians[3];
  unsigned int  size_structures[3], size_jacobians[3];
  unsigned int  start_structures[3], count_structures[3];
  unsigned long start[3], count[3];
  double        volumes[256], voxel_separations[3];
  int           i;
  double        *jacobians;
  double        *structures;
  double        volume;

  start[0] = 0;
  start[1] = 0;
  start[2] = 0;

  count[0] = 10;
  count[1] = 10;
  count[2] = 10;

  if (argc != 3) {
    fprintf(stderr, "Usage: label_volumes_from_jacobians structures.mnc jacobians.mnc\n");
    return(1);
  }

  /* read the structure volume */
  //printf("Opening %s\n", argv[1]);
  if (miopen_volume(argv[1], MI2_OPEN_READ, &hvol_structures) != MI_NOERROR) {
    fprintf(stderr, "Error opening volume: %s\n", argv[1]);
    return(1);
  }

  //printf("Opening %s\n", argv[2]);
  /* read the volume of jacobians */
  if (miopen_volume(argv[2], MI2_OPEN_READ, &hvol_jacobians) != MI_NOERROR) {
    fprintf(stderr, "Error opening volume: %s\n", argv[2]);
    return(1);
  }
  
  //printf("Getting Volume Dimensions\n");
  /* get dimensions */
  miget_volume_dimensions(hvol_structures, MI_DIMCLASS_SPATIAL, MI_DIMATTR_ALL,
			  MI_DIMORDER_FILE, 3, dimensions_structures);
  miget_volume_dimensions(hvol_jacobians, MI_DIMCLASS_SPATIAL, MI_DIMATTR_ALL,
			  MI_DIMORDER_FILE, 3, dimensions_jacobians);
  miget_dimension_sizes(dimensions_structures, 3, size_structures);
  miget_dimension_sizes(dimensions_jacobians, 3, size_jacobians);

  miget_dimension_separations(dimensions_structures, MI_ORDER_FILE, 
			      3, voxel_separations);


  for (i=0; i < 3; i++) {
    if (size_structures[i] != size_jacobians[i]) {
      fprintf(stderr, "Dimensions of two volumes must be the same.\n");
    }

    //printf("Width: %f\n", voxel_separations[i]);
  }

  count[0] = size_structures[0];
  count[1] = size_structures[1];
  count[2] = size_structures[2];

  //printf("Allocating memory\n");
  jacobians = (double *) malloc(size_jacobians[0] *
				size_jacobians[1] *
				size_jacobians[2] * sizeof(double));

  structures = (double *) malloc(size_structures[0] * 
			      size_structures[1] *
			      size_structures[2] * sizeof(double));


  //printf("Getting hyperslab.\n");
  if (miget_real_value_hyperslab(hvol_jacobians, MI_TYPE_DOUBLE,
				 (unsigned long *)start, 
				 (unsigned long *)count,
				 jacobians) < 0) {
    fprintf(stderr, "Could not get hyperslab\n");
  }
				 
  if (miget_real_value_hyperslab(hvol_structures, MI_TYPE_DOUBLE,
				 start, count,
				 structures) < 0) {
    fprintf(stderr, "Could not get hyperslab\n");
  }

  /* set all structure volumes to 0 */
  for (i=0; i < 255; i++) {
    volumes[i] = 0;
  }
  
  volume = voxel_separations[0] * voxel_separations[1] * voxel_separations[2];

  //printf("Computing tissue volumes.\n");
  for (i=0; i < count[0] * count[1] * count[2]; i++) {
    //printf("%i\n", (int)(structures[i] + 0.5));
    volumes[(int)(structures[i] + 0.5)] += volume * exp(jacobians[i]);
      
  }
  
  for (i=0; i < 255; i++) {
    if (volumes[i] != 0) {
      printf("%i, %f\n", i, volumes[i]);
    }
  }
    

  return(0);
}
