#!/usr/bin/env perl

use strict;
use warnings;
use MNI::Startup;
use MNI::Spawn;
use File::Basename;

my $usage = "$0 source_atlas target_atlas\n";

my $outputdir = undef;

my $source_atlas = shift or die $usage;
my $target_atlas = shift or die $usage;
$outputdir = shift;

if(defined($outputdir) and ! -d $outputdir)
{
	system("mkdir -p $outputdir");
}

my @source_blurs;
my @target_blurs;

my $source_base = "source";
my $target_base = "target";

if(!defined($outputdir))
{
	$outputdir = " ";
}

my $lsq12_1_xfm = "lsq12.xfm";
my $lsq12_2_xfm = "lsq12_2.xfm";
my $nlin_1_xfm = "nlin1.xfm";
my $nlin_2_xfm = "nlin2.xfm";
my $nlin_3_xfm = "nlin3.xfm";
my $nlin_4_xfm = "nlin4.xfm";

my @blurs = (0.25, 0.1);

foreach (@blurs) 
{
	Spawn("mincblur -clobber -gradient -fwhm $_ $source_atlas ${outputdir}${source_base}_$_");
	Spawn("mincblur -clobber -gradient -fwhm $_ $target_atlas ${outputdir}${target_base}_$_");
}

my $minctracc_opts = "-xcorr -clobber -debug ";

Spawn("minctracc $minctracc_opts -w_translations 0.2 0.2 0.2 -step 1 1 1 -lattice_diameter 3 3 3 -lsq12 -simplex 1 ${outputdir}${source_base}_0.25_blur.mnc ${outputdir}${target_base}_0.25_blur.mnc ${outputdir}${lsq12_1_xfm}");
Spawn("minctracc $minctracc_opts -w_translations 0.2 0.2 0.2 -step 1 1 1 -lattice_diameter 3 3 3 -lsq12 -simplex 1 -transform ${outputdir}lsq12.xfm ${outputdir}${source_base}_0.25_dxyz.mnc ${outputdir}${target_base}_0.25_dxyz.mnc ${outputdir}${lsq12_2_xfm}");

$minctracc_opts .= " -nonlinear corrcoeff -similarity 0.3 -stiffness 1 -weight 1 -sub_lattice 6 -w_translations 0.2 0.2 0.2";

Spawn("minctracc $minctracc_opts -step 0.5 0.5 0.5 -simplex 3 -iterations 60 -lattice_diameter 1.5 1.5 1.5 -transform ${outputdir}lsq12_2.xfm  ${outputdir}${source_base}_0.25_blur.mnc ${outputdir}${target_base}_0.25_blur.mnc ${outputdir}${nlin_1_xfm}");
Spawn("minctracc $minctracc_opts -step 0.5 0.5 0.5 -simplex 3 -iterations 60 -lattice_diameter 1.5 1.5 1.5 -transform ${outputdir}nlin1.xfm ${outputdir}${source_base}_0.25_dxyz.mnc ${outputdir}${target_base}_0.25_dxyz.mnc ${outputdir}${nlin_2_xfm}");
Spawn("minctracc $minctracc_opts -step 0.2 0.2 0.2 -simplex 1.5 -iterations 10 -lattice_diameter 0.6 0.6 0.6 -transform ${outputdir}nlin2.xfm ${outputdir}${source_base}_0.25_blur.mnc ${outputdir}${target_base}_0.25_blur.mnc ${outputdir}${nlin_3_xfm}");
Spawn("minctracc $minctracc_opts -step 0.2 0.2 0.2 -simplex 1.5 -iterations 10 -lattice_diameter 0.6 0.6 0.6 -transform ${outputdir}nlin3.xfm ${outputdir}${source_base}_0.25_dxyz.mnc ${outputdir}${target_base}_0.25_dxyz.mnc ${outputdir}${nlin_4_xfm}");

# Spawn("minctracc $minctracc_opts -step 0.1 0.1 0.1 -simplex 1 -iterations 10 -lattice_diameter 0.3 0.3 0.3 -transform nlin4.xfm $source_atlas $target_atlas nlin5.xfm");
