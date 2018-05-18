/* takes the dot product between a length of 3 deformation vector in a
   4D volume and the surface normals of a polyhedral object. 

   Author: Jason Lerch <jason@sickkids.ca>
*/
#define HAVE_MINC2 1

#include <volume_io.h>
#include <bicpl.h>
#include <stdio.h>

int main(int argc, char *argv[]) {
  VIO_Volume          grid_volume;
  VIO_File_formats    format;
  VIO_Real            dot_product, interp_values[3];
  object_struct   **objects;
  int             n_objects, sizes[VIO_MAX_DIMENSIONS], n_points, i;
  VIO_Vector          *normals;
  polygons_struct *polygons;
  VIO_STR          obj_filename, mnc_filename, output_filename;
  FILE            *file;
  VIO_Point           *points;

  /* get input arguments */
  initialize_argument_processing(argc, argv);
  if (!get_string_argument(NULL, &obj_filename) ||
      !get_string_argument(NULL, &mnc_filename) || 
      !get_string_argument(NULL, &output_filename)) {
    fprintf(stderr, "Usage: object_volume_dot_product polyhedra.obj displacement.mnc output.txt\n");
    return(1);
  }

  /* open the obj file */ 
  if (input_graphics_file(obj_filename,
			  &format, &n_objects, &objects ) != VIO_OK ) {
        return( 1 );
  }

  /* make sure data is polygons */
  if( n_objects != 1 || get_object_type(objects[0]) != POLYGONS ) {
        fprintf(stderr, "File must contain exactly 1 polygons struct.\n" );
        return( 1 );
  }

  polygons = get_polygons_ptr(objects[0]);
  n_points = get_object_points(objects[0], &points);
  normals = polygons->normals;

  /* open the displacement volume */
  if (input_volume(mnc_filename, 4, NULL, MI_ORIGINAL_TYPE, FALSE, 0.0, 0.0, 
		   TRUE, &grid_volume, (minc_input_options *)NULL) != VIO_OK) {
    return(1);
  }
  get_volume_sizes(grid_volume, sizes);

  /* open the output file */
  if( open_file(output_filename, WRITE_FILE, ASCII_FORMAT, &file ) != VIO_OK ) {
    return( 1 );
  }

  for (i = 0; i < n_points; i++) { 
    evaluate_volume_in_world(grid_volume, Point_x(points[i]), 
			     Point_y(points[i]), Point_z(points[i]),
			     -1, TRUE, 0.0, interp_values,
			     NULL, NULL, NULL, NULL, NULL, NULL,
			     NULL, NULL, NULL);

    dot_product =  Vector_x(normals[i]) * interp_values[0];
    dot_product += Vector_y(normals[i]) * interp_values[1];
    dot_product += Vector_z(normals[i]) * interp_values[2];


    /*
    printf("%f %f %f %f %f %f %f\n", Vector_x(normals[i]),
	   Vector_y(normals[i]), Vector_z(normals[i]),
	   interp_values[0], interp_values[1], interp_values[2],
	   dot_product);
    */

    (void) output_real(file, dot_product);
    (void) output_newline(file);
  }
  
  (void) close_file(file);
  return(0);
}
    
