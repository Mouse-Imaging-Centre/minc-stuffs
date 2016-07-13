/* ----------------------------- MNI Header -----------------------------------
@NAME       : tagtoxfm
@INPUT      : argc, argv - command line arguments
@OUTPUT     : (none)
@RETURNS    : status
@DESCRIPTION: Program to calculate a transform file from a tag file.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : August 30, 1993 (Peter Neelin)
@MODIFIED   : $Log: tagtoxfm.c,v $
@MODIFIED   : Revision 1.7  2001-05-23 04:13:04  stever
@MODIFIED   : Merge from branch-1_3 branch.
@MODIFIED   :
@MODIFIED   : Revision 1.6.2.1  2000/11/12 15:49:54  stever
@MODIFIED   : Removing obsolete Build directory
@MODIFIED   :
@MODIFIED   : Revision 1.6  1999/06/21 20:18:19  stever
@MODIFIED   : final checkin before switch to CVS
@MODIFIED   :
 * Revision 1.5  1997/12/10  20:25:06  david
 * check_in_all
 *
 * Revision 1.4  1995/12/19  15:47:00  david
 * check_in_all
 *
 * Revision 1.3  1995/12/18  16:43:40  david
 * *** empty log message ***
 *
 * Revision 1.7  1995/12/15  21:29:25  neelin
 * Recompiled (making modifications for changes to volume_io) so that
 * calculations are done in double precision.
 *
 * Revision 1.6  94/04/22  08:17:26  neelin
 * Changed back to using compute_transform_from_tags.
 * 
 * Revision 1.5  94/04/22  08:08:40  neelin
 * Modified to use safe_compute_transform_from_tags and volume_io.h instead
 * of def_mni.h
 * 
 * Revision 1.4  93/09/16  10:02:45  neelin
 * Added open_file_with_default_suffix for reading tag files and 
 * output_transform_file for writing xfm files to support default suffixes.
 * 
 * Revision 1.3  93/09/08  15:45:42  neelin
 * Added checking for number of points and writing of comments in .xfm file.
 * 
 * Revision 1.2  93/09/01  15:36:15  neelin
 * Changed usage error message.
 * 
 * Revision 1.1  93/09/01  15:25:33  neelin
 * Initial revision
 * 
 * Revision 1.1  93/09/01  13:22:02  neelin
 * Initial revision
 * 
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /static-cvsroot/visualization/Register/Tagtoxfm/tagtoxfm.c,v 1.7 2001-05-23 04:13:04 stever Exp $";
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* in order to compile against the minc-toolkit without CMake */
#define HAVE_MINC2 1
#include <volume_io.h>
/* \in order to... */

#include <bicpl.h>
#include <ParseArgv.h>
#include <bicpl/compute_xfm.h>
#include "tagtoxfm.h"

/* Constants */
#ifndef TRUE
#  define TRUE 1
#  define FALSE 0
#endif

/* Main program */

int main(int argc, char *argv[])
{
   char *pname, *tagfile, *xfmfile;
   int n_volumes, n_tag_points;
   VIO_Real **tags_volume1, **tags_volume2, **temp_tags;
   VIO_General_transform transform;
   FILE *fp;
   char *type_string, *inverse_string, comment[512];

   /* Parse arguments */
   pname = argv[0];
   if (ParseArgv(&argc, argv, argTable, 0) || (argc != 3)) {
      (void) fprintf(stderr, 
                     "\nUsage: %s [<options>] infile.tag outfile.xfm\n\n",
                     argv[0]);
      exit(EXIT_FAILURE);
   }
   tagfile = argv[1];
   xfmfile = argv[2];

   /* Read in tag file */
   if ((open_file_with_default_suffix(tagfile,
                  get_default_tag_file_suffix(),
                  READ_FILE, ASCII_FORMAT, &fp) != VIO_OK) ||
       (input_tag_points(fp, &n_volumes, &n_tag_points, 
                         &tags_volume1, &tags_volume2, 
                         NULL, NULL, NULL, NULL) != VIO_OK)) {
      (void) fprintf(stderr, "%s: Error reading tag file %s\n", 
                     pname, tagfile);
      exit(EXIT_FAILURE);
   }
   (void) close_file(fp);

   /* Check number of volumes */
   if (n_volumes != 2) {
      (void) fprintf(stderr, "%s: Wrong number of volumes in %s\n", 
                     pname, tagfile);
      exit(EXIT_FAILURE);
   }

   /* Check number of points for linear transformation */
   if (((transform_type == TRANS_LSQ6) ||
        (transform_type == TRANS_LSQ7) ||
        (transform_type == TRANS_LSQ9) ||
        (transform_type == TRANS_LSQ10) ||
        (transform_type == TRANS_LSQ12)) &&
       (n_tag_points < MIN_POINTS_LINEAR) ) {
      (void) fprintf(stderr, 
                     "%s: Need at least %d points (only %d in %s)\n", 
                     pname, MIN_POINTS_LINEAR, n_tag_points, tagfile);
      exit(EXIT_FAILURE);
   }

   /* Check number of points for thin-plate spline transformation */
   if ((transform_type == TRANS_TPS) &&
       (n_tag_points < MIN_POINTS_TPS) ) {
      (void) fprintf(stderr, 
                     "%s: Need at least %d points (only %d in %s)\n", 
                     pname, MIN_POINTS_TPS, n_tag_points, tagfile);
      exit(EXIT_FAILURE);
   }

   /* If inverting, switch order of points */
   if (inverse) {
      temp_tags = tags_volume1;
      tags_volume1 = tags_volume2;
      tags_volume2 = temp_tags;
   }

   /* Compute transformation */
   compute_transform_from_tags(n_tag_points, tags_volume1, tags_volume2,
                               transform_type, &transform);

   /* Create output file comment */
   switch (transform_type) {
   case TRANS_LSQ6: type_string = "6 parameter linear least-squares"; break;
   case TRANS_LSQ7: type_string = "7 parameter linear least-squares"; break;
   case TRANS_LSQ9: type_string = "9 parameter linear least-squares"; break;
   case TRANS_LSQ10: type_string = "10 parameter linear least-squares"; break;
   case TRANS_LSQ12: type_string = "12 parameter linear least-squares"; break;
   case TRANS_TPS: type_string = "thin-plate spline"; break;
   default: type_string = "unknown"; break;
   }
   if (inverse)
      inverse_string = " with -inverse";
   else
      inverse_string = "";
   (void) sprintf(comment, " Created from tag file %s\n using %s%s",
                  tagfile, type_string, inverse_string);

   /* Save transformation */
   if (!clobber && file_exists(xfmfile)) {
      (void) fprintf(stderr, "%s: xfm file \"%s\" already exists.\n",
                     pname, xfmfile);
      exit(EXIT_FAILURE);
   }
   if (output_transform_file(xfmfile, comment, &transform) != VIO_OK) {
      (void) fprintf(stderr, "%s: Error writing xfm file %s\n", 
                     pname, xfmfile);
      exit(EXIT_FAILURE);
   }

   return EXIT_SUCCESS;
}

