#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <volume_io.h>
#include <bicpl.h>
#include <ParseArgv.h>
#include <time_stamp.h>

/* argument parsing defaults */
static int verbose = FALSE;
static int clobber = FALSE;

/* argument table */
static ArgvInfo argTable[] = {
   {NULL, ARGV_HELP, (char *)NULL, (char *)NULL,
    "General options:"},
   {"-verbose", ARGV_CONSTANT, (char *)TRUE, (char *)&verbose,
    "print out extra information"},
   {"-clobber", ARGV_CONSTANT, (char *)TRUE, (char *)&clobber,
    "clobber existing files"},
   {NULL, ARGV_HELP, (char *)NULL, (char *)NULL, ""},
   {NULL, ARGV_END, NULL, NULL, NULL}
};

int main(int argc, char *argv[]) {
  int v1, v2, v3, v4;
  int sizes[MAX_DIMENSIONS], grid_sizes[4];
  int n_concat_transforms, i;
  char *arg_string;
  char *input_volume_name;
  char *input_xfm;
  char *outfile;
  Real w1, w2, w3;
  Real nw1, nw2, nw3;
  Real original[3], transformed[3];
  Real value;
  Real cosine[3];
  Real original_separation[3], grid_separation[4];
  Real original_starts[3], grid_starts[4];
  Volume eval_volume, new_grid;
  General_transform xfm, *voxel_to_world;
  STRING *dimnames, dimnames_grid[4];
  progress_struct progress;
  
  arg_string = time_stamp(argc, argv);

  /* Check arguments   */
  if(ParseArgv(&argc, argv, argTable, 0) || (argc != 4)){
    fprintf(stderr, 
	    "\nUsage: %s [options] input.mnc input.xfm output_grid.mnc\n", 
	    argv[0]);
    fprintf(stderr, "       %s -help\n\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  input_volume_name = argv[1];
  input_xfm = argv[2];
  outfile = argv[3];

  /* check for the infile and outfile */
  if(access(input_volume_name, F_OK) != 0){
    fprintf(stderr, "%s: Couldn't find %s\n\n", argv[0], input_volume_name);
    exit(EXIT_FAILURE);
  }
  if(access(input_xfm, F_OK) != 0) {
    fprintf(stderr, "%s: Couldn't find %s\n\n", argv[0], input_xfm);
    exit(EXIT_FAILURE);
  }
  if(access(outfile, F_OK) == 0 && !clobber){
    fprintf(stderr, "%s: %s exists! (use -clobber to overwrite)\n\n", 
	    argv[0], outfile);
    exit(EXIT_FAILURE);
  }

  /*--- input the volume */
  /*
  if( input_volume( input_volume_name, 3, NULL, MI_ORIGINAL_TYPE, 
                    FALSE, 0.0, 0.0, TRUE, &eval_volume,
		    (minc_input_options *) NULL ) != OK )
    return( 1 );
  */

  if (input_volume_header_only( input_volume_name, 3, NULL, &eval_volume,
				(minc_input_options *) NULL ) != OK ) 
    return( 1 );

  /* get information about the volume */
  get_volume_sizes( eval_volume, sizes );
  voxel_to_world = get_voxel_to_world_transform(eval_volume);
  dimnames = get_volume_dimension_names(eval_volume);
  get_volume_separations(eval_volume, original_separation);
  get_volume_starts(eval_volume, original_starts);

  /* create new 4D volume, last three dims same as other volume,
     first dimension being the vector dimension. */
  for(i=1; i < 4; i++) {
    dimnames_grid[i] = dimnames[i-1];
    grid_separation[i] = original_separation[i-1];
    grid_sizes[i] = sizes[i-1];
    grid_starts[i] = original_starts[i-1];
  }
  dimnames_grid[0] = "vector_dimension";
  grid_sizes[0] = 3;
  grid_separation[0] = 1;
  grid_starts[0] = 0;

  new_grid = create_volume(4, dimnames_grid, NC_SHORT, FALSE, 0.0, 0.0);

  //set_voxel_to_world_transform(new_grid, voxel_to_world);
  // initialize the new grid volume, otherwise the output will be
  // garbage...
  set_volume_real_range(new_grid, -100, 100);
  set_volume_sizes(new_grid, grid_sizes);

  set_volume_separations(new_grid, grid_separation);
  set_volume_starts(new_grid, grid_starts);
  /*
  for (i=0; i < 3; i++) {
    get_volume_direction_cosine(eval_volume, i, cosine);
    set_volume_direction_cosine(new_grid, i+1, cosine);
  }
  */

  alloc_volume_data(new_grid);

  /* get the transforms */
  if( input_transform_file( input_xfm, &xfm ) != OK )
    return( 1 );

  /* see how many transforms will be applied */
  n_concat_transforms = get_n_concated_transforms( &xfm );
  printf("Number of transforms to be applied: %d\n", n_concat_transforms);

  initialize_progress_report(&progress, FALSE, sizes[0], "Processing");


  /* evaluate the transform at every voxel, keep the displacement
     in the three cardinal directions */
  for( v1 = 0;  v1 < sizes[0];  ++v1 ) {
    update_progress_report(&progress, v1 + 1);
    for( v2 = 0;  v2 < sizes[1];  ++v2 ) {
      for( v3 = 0;  v3 < sizes[2];  ++v3 ) {
	convert_3D_voxel_to_world(eval_volume, v1, v2, v3, &original[0], 
				  &original[1], &original[2]);
	general_transform_point(&xfm, original[0], original[1], original[2], 
				&transformed[0], &transformed[1], 
				&transformed[2]);
	for(i=0; i < 3; i++) {
	  value = transformed[i] - original[i];
	  set_volume_real_value(new_grid, i, v1, v2, v3, 0, value);
	}
      }
    }
  }

  terminate_progress_report(&progress);

  printf("Outputting volume.\n");

  output_volume(outfile, MI_ORIGINAL_TYPE, TRUE, 0.0, 0.0, 
		new_grid, arg_string, NULL);

  return(0);

}
