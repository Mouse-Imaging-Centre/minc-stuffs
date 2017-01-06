#!/usr/bin/perl

# generates an xfm file for a provided grid file

use strict;
use warnings;

use File::Spec;
use MNI::Startup;
use Getopt::Tabular;
use MNI::Spawn;
use MNI::FileUtilities qw(test_file check_output_dirs);

my $usage = "

$0 [options] deformation_grid.mnc output_xfm_for_grid.xfm\n";


# handle arguments
my @left_over_args;
my @arg_table = 
    ( @DefaultArgs,
    );

GetOptions(\@arg_table, \@ARGV, \@left_over_args) or die "\n";


my $deformation_grid = shift @left_over_args or die $usage;
my $output_xfm = shift @left_over_args or die $usage;

if(-f $output_xfm and not $Clobber) {
  die "Error: output xfm file already exists and -clobber not specified: $output_xfm.\n";
}

if(not -f $deformation_grid){
  die "Error: deformation grid does not exist: $deformation_grid.\n";
}

# save the realpath to the grid, otherwise it won't necessarily work
my $deformation_grid_abs_path = File::Spec->rel2abs( $deformation_grid );

# write to the output
open(OUTPUT, "> $output_xfm");
print OUTPUT "MNI Transform File\n";
print OUTPUT "% created by $0 with grid $deformation_grid\n";
print OUTPUT "Transform_Type = Grid_Transform;\n";
print OUTPUT "Displacement_Volume = $deformation_grid_abs_path;\n";
