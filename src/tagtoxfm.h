/* ----------------------------- MNI Header -----------------------------------
@NAME       : tagtoxfm.h
@DESCRIPTION: Header file for tagtoxfm
@GLOBALS    : 
@CALLS      : 
@CREATED    : August 30, 1993 (Peter Neelin, Louis Collins)
@MODIFIED   : 
---------------------------------------------------------------------------- */

/* Argument variables */

/*  These are not needed?  - DMD, Aug. 25, 1998
extern double ftol;
extern double simplex_size;
*/

int 
  inverse = FALSE;
Trans_type
  transform_type = TRANS_LSQ6;
int clobber = FALSE;

/* Argument table */
ArgvInfo argTable[] = {
   {NULL, ARGV_HELP, NULL, NULL,
       "---Transformation maps volume two to volume one---"},
   {NULL, ARGV_HELP, NULL, NULL,
       "Transformation type. Default = -lsq6."},
   {"-lsq6", ARGV_CONSTANT, (char *) TRANS_LSQ6, (char *) &transform_type,
       "6 parameter (scale=1.0) least-squares linear transformation."},
   {"-lsq7", ARGV_CONSTANT, (char *) TRANS_LSQ7, (char *) &transform_type,
       "7 parameter (one scale) least-squares linear transformation."},
   {"-lsq9", ARGV_CONSTANT, (char *) TRANS_LSQ9, (char *) &transform_type,
       "9 parameter least-squares linear transformation."},
   {"-lsq10", ARGV_CONSTANT, (char *) TRANS_LSQ10, (char *) &transform_type,
       "10 parameter least-squares linear transformation."},
   {"-lsq12", ARGV_CONSTANT, (char *) TRANS_LSQ12, (char *) &transform_type,
       "12 parameter least-squares linear transformation."},
   {"-tps", ARGV_CONSTANT, (char *) TRANS_TPS, (char *) &transform_type,
       "Thin-plate spline non-linear transformation."},
   {NULL, ARGV_HELP, NULL, NULL,
       "Other options:"},
   {"-inverse", ARGV_CONSTANT, (char *) TRUE, (char *) &inverse,
       "Swap tags, then compute transform (default=FALSE)."},
#ifdef NOT_NEEDED
   {"-tol", ARGV_FLOAT, (char *) 0,(char *) &ftol,
       "Stopping criteria tolerence"},
   {"-simplex", ARGV_FLOAT, (char *) 0, (char *) &simplex_size,
       "Radius of simplex volume."},
#endif
   {"-clobber", ARGV_CONSTANT, (char *) TRUE, (char *) &clobber,
       "Overwrite any existing xfm file."},
   {"-noclobber", ARGV_CONSTANT, (char *) FALSE, (char *) &clobber,
       "Do not overwrite any existing xfm file."},

   {NULL, ARGV_END, NULL, NULL, NULL}
};

