/* Takes a transform and a volume as the input, and outputs a tag file
   which can be used by tagtoxfm. The primary goal is to compute a
   linear transformation from a non-linear grid transform.
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define HAVE_MINC2 1

#include <volume_io.h>
#include <bicpl.h>
#include <ParseArgv.h>
#include <time_stamp.h>

/* argument parsing defaults */
static int verbose = FALSE;
static int clobber = FALSE;

char *mask_file;

/* argument table */
static ArgvInfo argTable[] = {
   {NULL, ARGV_HELP, (char *)NULL, (char *)NULL,
    "General options:"},
   {"-verbose", ARGV_CONSTANT, (char *)TRUE, (char *)&verbose,
    "print out extra information"},
   {"-clobber", ARGV_CONSTANT, (char *)TRUE, (char *)&clobber,
    "clobber existing files"},
   {"-mask", ARGV_STRING, (char *)1, (char *)&mask_file,
    "<mask.mnc> Use mask file for calculations."},

   {NULL, ARGV_HELP, (char *)NULL, (char *)NULL, ""},
   {NULL, ARGV_END, NULL, NULL, NULL}
};

int main(int argc, char *argv[]) {
  VIO_Volume            eval_volume, mask_volume;
  VIO_General_transform xfm;
  char                  *output_tagfile_name;
  char                  *input_volume_name;
  char                  *input_transform_name;
  int                   sizes[VIO_MAX_DIMENSIONS];
  FILE                  *tagfile;
  VIO_Real              world_coordinates[3];
  VIO_Real              target_coordinates[3];
  int                   counter_offset, use_voxel;
  int                   v1,v2,v3;

  /* handle arguments */
  if(ParseArgv(&argc, argv, argTable, 0) || (argc != 4)){
    fprintf(stderr,
            "\nUsage: %s [options] input.mnc input.xfm output_tags.tag\n",
            argv[0]);
    fprintf(stderr, "       %s -help\n\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  /* get filenames */
  input_volume_name    = argv[1];
  input_transform_name = argv[2];
  output_tagfile_name  = argv[3];

  /* check for the infile and outfile */
  if(access(input_volume_name, F_OK) != 0){
    fprintf(stderr, "%s: Couldn't find %s\n\n", argv[0], input_volume_name);
    exit(EXIT_FAILURE);
  }
  if(access(input_transform_name, F_OK) != 0) {
    fprintf(stderr, "%s: Couldn't find %s\n\n", argv[0], input_transform_name);
    exit(EXIT_FAILURE);
  }
  if(access(output_tagfile_name, F_OK) == 0 && !clobber){
    fprintf(stderr, "%s: %s exists! (use -clobber to overwrite)\n\n",
            argv[0], output_tagfile_name);
    exit(EXIT_FAILURE);
  }

  /* input the volume */
  if( input_volume_header_only( input_volume_name, 3, NULL, &eval_volume,
				(minc_input_options *) NULL ) != VIO_OK )
    return( 1 );

  /* input the mask if so desired*/
  if (mask_file != NULL) {
    if( input_volume( mask_file, 3, NULL, MI_ORIGINAL_TYPE,
		      FALSE, 0.0, 0.0, TRUE, &mask_volume, 
		      (minc_input_options *) NULL) != VIO_OK )
      return( 1 );
  }

  /* get information about the volume */
  get_volume_sizes( eval_volume, sizes );
  
  /* get the transforms */
  if( input_transform_file( input_transform_name, &xfm ) != VIO_OK )
    return( 1 );

  /* open the tagfile for editing */
  tagfile = fopen(output_tagfile_name, "w");
  initialize_tag_file_output( tagfile, "Tagfile from xfm", 2 );
  
  counter_offset = 5;
  for (v1=0; v1 < (sizes[0] - counter_offset); v1 += counter_offset) {
    for (v2=0; v2 < (sizes[1] - counter_offset); v2 += counter_offset) {
      for (v3=0; v3 < (sizes[2] - counter_offset); v3 += counter_offset) {
	use_voxel = 1;

	/* if mask file exists check whether current voxel is inside mask */
	if (mask_file != NULL) {
	  if (get_volume_real_value(mask_volume, v1, v2, v3, 0, 0) < 0.5) {
	    use_voxel = 0;
	  }
	}

	if (use_voxel == 1) { /* will always be true in absence of mask */
	    
	  convert_3D_voxel_to_world(eval_volume, v1, v2, v3, 
				    &world_coordinates[0], 
				    &world_coordinates[1],
				    &world_coordinates[2]);
	  general_transform_point(&xfm, world_coordinates[0],
				  world_coordinates[1], world_coordinates[2],
				  &target_coordinates[0], 
				  &target_coordinates[1],
				  &target_coordinates[2]);
	  
	  output_one_tag(tagfile, 2, target_coordinates, world_coordinates,
			 NULL, NULL, NULL, NULL);
	}
      }
    }
  }
  terminate_tag_file_output(tagfile);
  fclose(tagfile);
  return(0);
}
