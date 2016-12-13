#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pyminc.volumes.factory import *
from numpy import *
from argparse import ArgumentParser
import tempfile
import os
import functools
import signal
import shutil
import subprocess
import sys
import math


def get_tempfile(suffix):
    #tmp_fd, tmpfile = tempfile.mkstemp(suffix='.xfm')
    #os.close(tmp_fd)
    counter = 0
    pid = os.getpid()
    tmpdir = "%s/rot_%s" % (os.environ["TMPDIR"], pid)
    if not os.access(tmpdir, 0):
        try:
            subprocess.check_call(("mkdir -p %s" % tmpdir).split())
        except:
            sys.exit("ERROR: could not make temporary directory %s!" % tmpdir)
            
    tmpfile = "/%s/rot_%s.%s" % (tmpdir, counter, suffix)
    while os.access(tmpfile, 0):
        counter = counter + 1
        tmpfile = "/%s/rot_%s%s" % (tmpdir, counter, suffix)
    return tmpfile

def get_centre_of_gravity(file):
    lines = subprocess.check_output(["volume_cog", file])
    line = lines.decode().split('\n')[-2]
    cog = array(line.strip().split(" ")).astype("float32")
    return cog

def get_coordinates_from_tag_file(tag_file):
    #
    # Tag file content looks as follows:
    #
    # MNI Tag Point File
    # Volumes = 2;
    # %Mon Apr 25 16:58:13 2016>>> find_peaks -pos_only -min_distance 2 //dev/shm//rot_23705/rot_5.mnc //dev/shm//rot_23705/rot_0..tag
    #
    # Points =
    # -28.6989 -49.7709 -28.283 -28.6989 -49.7709 -28.283 "12"
    # -28.5549 -49.7649 -25.445 -28.5549 -49.7649 -25.445 "11"
    # -28.1424 -49.6824 -13.28 -28.1424 -49.6824 -13.28 "11"
    # -27.7989 -49.2069 -22.667 -27.7989 -49.2069 -22.667 "10";
    #
    all_coordinates = []
    with open(tag_file) as f:
        content = f.readlines()
        found_coordinates = False
        for line in content:
            if not found_coordinates:
                if "Points" in line:
                    found_coordinates = True
            else:
                # each line with coordinates starts with a space
                # so the very first element when splitting this 
                # line is the empty string. We don't need that one
                current_coor = array([line.split(' ')[1],
                                      line.split(' ')[2],
                                      line.split(' ')[3]]).astype("float32")
                all_coordinates.append(current_coor)
    return all_coordinates
    
    
def get_distance_transform_peaks(input_file, peak_distance):
    #
    # for solid input files (brains, embryos) we can
    # calculate the distance transform for the input file
    # and use peaks from that distance transform to 
    # widen our search space a bit.
    #
    # the .decode() in the end removes the b from the front
    # of the returned string (Python3 default)
    bimodalt_value = subprocess.check_output(["mincstats",
                                             "-quiet",
                                             "-biModalT",
                                             input_file]).rstrip().decode()
    max_value = subprocess.check_output(["mincstats",
                                         "-quiet",
                                         "-max",
                                         input_file]).rstrip().decode()
    distance_transform = get_tempfile('.mnc')
    subprocess.check_call(("mincmorph -successive B[%s:%s]F %s %s" %
                            (bimodalt_value,max_value,input_file,distance_transform)).split())
    peak_tags = get_tempfile('.tag')
    subprocess.check_call(("find_peaks -pos_only -min_distance %s %s %s" % (peak_distance, distance_transform, peak_tags)).split())
    all_coors = get_coordinates_from_tag_file(peak_tags)
    return all_coors
    
def get_blur_peaks(input_file, blur_kernel, peak_distance):
    blurred_input = get_tempfile('_blur.mnc')
    subprocess.check_call(("mincblur -no_apo -fwhm %s %s %s" % (blur_kernel, input_file, blurred_input.split('_blur.mnc')[0])).split())
    peak_tags = get_tempfile('.tag')
    subprocess.check_call(("find_peaks -pos_only -min_distance %s %s %s" % (peak_distance, blurred_input, peak_tags)).split())
    all_coors = get_coordinates_from_tag_file(peak_tags)
    return all_coors

def compute_xcorr(sourcefile, targetvol, maskvol):
    try:
        sourcevol = volumeFromFile(sourcefile)
    except mincException:
        raise
    #    return 0
    
    f1 = 0
    f2 = 0
    f3 = 0
    
    if maskvol is not None:
        # we've had issues with mask files that after 
        # running autocrop on them, only contain "nan"
        # values. This seems to happen when the original
        # mask has a valid image range indicated as:
        # image: unsigned short 1 to 1
        # when that happens, all following calculations 
        # fail, so we should quit:
        if math.isnan(maskvol.data.sum()):
            # clean up...
            shutil.rmtree("%s/rot_%s" % (os.environ["TMPDIR"], os.getpid()))
            raise ValueError("\n\n* * * * * * * * * *\n"
                  "Error: the mask volume is corrupted. No values are found inside the mask.\n"
                  "probably you are using a mask with only the value 1 in it, but which at the same\n"
                  "time has a valid range of 1 to 1 . You can check this using mincinfo.\n"
                  "You need to generate a new mask. To produce a full field of view mask for file.mnc\n"
                  "run:\n\nminccalc -expression if(true){out = 1;} file.mnc file_mask.mnc\n"
                  "* * * * * * * * * *\n")
        f1 = sum(sourcevol.data[maskvol.data > 0.5] * \
             targetvol.data[maskvol.data > 0.5])
        f2 = sum(sourcevol.data[maskvol.data > 0.5] ** 2)
        f3 = sum(targetvol.data[maskvol.data > 0.5] ** 2)
    else:
        f1 = sum(sourcevol.data * targetvol.data)
        f2 = sum(sourcevol.data ** 2)
        f3 = sum(targetvol.data ** 2)

    sourcevol.closeVolume()

    return f1 / (sqrt(f2) * sqrt(f3))

def create_transform(cog_diff, xrot, yrot, zrot, cog_source):
    # create temporary file for transform
    tmp_transform = get_tempfile('.xfm')
    subprocess.check_call(("param2xfm -translation %s %s %s -rotation %s %s %s -center %s %s %s %s"
                           % (cog_diff[0], cog_diff[1], cog_diff[2],
                              xrot, yrot, zrot, 
                              cog_source[0], cog_source[1], cog_source[2],
                              tmp_transform)).split())
    return tmp_transform

def resample_volume(source, target, transform):
    tmp_resampled = get_tempfile('.mnc')
    subprocess.check_call(("mincresample -transform %s -like %s %s %s" 
                          % (transform, target, source, tmp_resampled)).split())
    return tmp_resampled

def minctracc(source, target, mask, stepsize, wtranslations, simplex):
    wtrans_decomp = array(wtranslations.split(',')).astype("float")
    tmp_transform = get_tempfile('.xfm')
    cmd = ("minctracc -identity -lsq6 -xcorr -simplex %s -step %s %s %s %s %s %s -w_translations %s %s %s " 
           % (simplex, stepsize, stepsize, stepsize, source, target, tmp_transform,
           wtrans_decomp[0], wtrans_decomp[1], wtrans_decomp[2]))
    if mask:
        cmd += ("-source_mask %s -model_mask %s " % (mask, mask))
    print(cmd)
    subprocess.check_call(cmd.split())
    
    return tmp_transform

def concat_transforms(t1, t2):
    tmp_transform = get_tempfile('.xfm')
    subprocess.check_call(("xfmconcat %s %s %s" % (t1, t2, tmp_transform)).split())
    return tmp_transform

def get_cross_correlation_from_coordinate_pair(source_img, target_img, target_vol, mask, coordinate_pair):
    # Generate a transformation based on the coordinate_pair provided. Apply this 
    # transformation to the source_img and calculate the cross correlation between
    # the source_img and target_img after this initial alignment
    #
    #
    # source coordinates : coordinate_pair[0]
    # target_coordinates : coordinate_pair[1]
    transform_from_coordinates = create_transform(coordinate_pair[1] - coordinate_pair[0], 0, 0, 0, coordinate_pair[0])
    resampled_source = resample_volume(source_img, target_img, transform_from_coordinates)
    xcorr = compute_xcorr(resampled_source, target_vol, mask)
    os.remove(resampled_source)
    return float(xcorr)
    
def loop_rotations(stepsize, source, target, mask, simplex, start=50, interval=10, 
                   wtranslations="0.2,0.2,0.2", use_multiple_seeds=True, max_number_seeds=5):
    # load the target and mask volumes
    targetvol = volumeFromFile(target)
    maskvol = volumeFromFile(mask) if mask is not None else None
    
    # 1) The default way of aligning files is by using the centre
    #    of gravity of the input files
    cog_source = get_centre_of_gravity(source)
    cog_target = get_centre_of_gravity(target)
    list_of_coordinate_pairs = [[cog_source, cog_target]]
    
    # 2) If we are using multiple seeds, calculate possible
    #    seeds for both source and target images. The distance
    #    between peaks is based on the stepsize used for the
    #    registrations
    if use_multiple_seeds:
        list_source_peaks = get_distance_transform_peaks(input_file=source, peak_distance=stepsize)
        print("\n\nPeaks found in the source image (Distance Transform):")
        for coor_src in list_source_peaks:
            print(coor_src)
        # also add peaks from the blurred version of the input file
        blurred_peaks_source = get_blur_peaks(input_file=source, blur_kernel=stepsize, peak_distance=stepsize)
        print("\n\nPeaks found in the source image (blurrred image):")
        for coor_src in blurred_peaks_source:
            print(coor_src)
            list_source_peaks.append(coor_src)
        # also add the center of gravity of the source image
        list_source_peaks.append(cog_source)
        list_target_peaks = get_distance_transform_peaks(input_file=target, peak_distance=stepsize)
        print("\n\nPeaks found in the target image (Distance Transform):")
        for coor_trgt in list_target_peaks:
            print(coor_trgt)
        blurred_peaks_target = get_blur_peaks(input_file=target, blur_kernel=stepsize, peak_distance=stepsize)
        print("\n\nPeaks found in the target image (blurrred image):")
        for coor_target in blurred_peaks_target:
            print(coor_target)
            list_target_peaks.append(coor_target)
        # same for the target; add the center of gravity:
        list_target_peaks.append(cog_target)
        for source_coor in list_source_peaks:
            for target_coor in list_target_peaks:
                list_of_coordinate_pairs.append([source_coor, target_coor])
    
        # 3) If we have more coordinates pairs than we'll be using, we'll have
        #    to determine the initial cross correlation of each pair and sort
        #    the pairs based on that
        pairs_with_xcorr = []
        if len(list_of_coordinate_pairs) > max_number_seeds:
            for coor_pair in list_of_coordinate_pairs:
                xcorr_coor_pair = get_cross_correlation_from_coordinate_pair(source, target, targetvol, maskvol, coor_pair)
                pairs_with_xcorr.append({'xcorr': xcorr_coor_pair,
                                         'coorpair': coor_pair})
                print("Xcorr: " + str(xcorr_coor_pair))
            sort_results(pairs_with_xcorr, reverse_order=True)
            print("\n\n\nCoordinate pairs and their cross correlation: \n\n")
            print(pairs_with_xcorr)
            # empty the list and fill it up with the best matches:
            list_of_coordinate_pairs = []
            for i in range(max_number_seeds):
                list_of_coordinate_pairs.append(pairs_with_xcorr[i]['coorpair'])
            print("\n\nNew list of coordinates:")
            print(list_of_coordinate_pairs)
    
    results = []
    best_xcorr = 0
    for coordinates_src_target in list_of_coordinate_pairs:
        coor_src = coordinates_src_target[0]
        coor_trgt = coordinates_src_target[1]
        for x in range(-start, start+1, interval):
            for y in range(-start, start+1, interval):
                for z in range(-start, start+1, interval):
                    # we need to include the centre of the volume as rotation centre = cog1
                    init_transform = create_transform(coor_trgt - coor_src, x, y, z, coor_src)
                    init_resampled = resample_volume(source, target, init_transform)
                    transform = minctracc(init_resampled, target, mask, stepsize=stepsize,
                                          wtranslations=wtranslations, simplex=simplex)
                    resampled = resample_volume(init_resampled, target, transform)
                    conc_transform = concat_transforms(init_transform, transform)
                    xcorr = compute_xcorr(resampled, targetvol, maskvol)
                    if isnan(xcorr):
                        xcorr = 0
                    results.append({'xcorr': xcorr, 'transform': conc_transform, \
                                    'resampled': resampled, 'x': x, \
                                    'y': y, 'z': z})
                    if xcorr > best_xcorr:
                        best_xcorr = xcorr
                    # had some issues with the resampled file being gone...
                    # we'll just resample the final file only at the end
                    os.remove(resampled)
                    os.remove(init_resampled)
                    print("FINISHED: %s %s %s :: %s" % (x,y,z, xcorr))
    
    sort_results(results)
    # resample the best result:
    final_resampled = resample_volume(source, target, results[-1]["transform"])
    results[-1]["resampled"] = final_resampled
    targetvol.closeVolume()
    if mask is not None:
        maskvol.closeVolume()
    return results

def dict_extract(adict, akey):
    return adict[akey]

def sort_results(results, reverse_order=False):
    sort_key_func = functools.partial(dict_extract, akey="xcorr")
    results.sort(key=sort_key_func, reverse=reverse_order)

def downsample(infile, stepsize):
    output = get_tempfile(".mnc")
    subprocess.check_call(("autocrop -isostep %s %s %s" % (stepsize, infile, output)).split())
    return output

# clean up tmp on soft kill signal
def termtrapper(signum, frame):
    print("got kill signal")
    os.removedirs("%s/rot_%s" % (os.environ["TMPDIR"], os.getpid()))
    exit(1)

if __name__ == "__main__":
    # handle soft kill signal to clean up tmp
    signal.signal(signal.SIGTERM, termtrapper)

    parser = ArgumentParser()
    parser.add_argument("-m", "--mask", dest="mask",
                      help="mask to use for computing xcorr", type=str)
    parser.add_argument("-s", "--stepsizeresample", dest="resamplestepsize",
                      help="resampled volumes to this stepsize [default = %(default)s]",
                      type=float, default=0.2)
    parser.add_argument("-g", "--stepsizeregistration", dest="registrationstepsize",
                      help="use this stepsize in the minctracc registration [default = %(default)s]",
                      type=float, default=0.6)
    parser.add_argument("-t", "--tempdir", dest="tmpdir",
                      help="temporary directory to use",
                      type=str)
    parser.add_argument("-r", "--range", dest="range",
                      help="range of rotations to search across [default = %(default)s]",
                      type=int, default=50)
    parser.add_argument("-i", "--interval", dest="interval",
                      help="interval (in degrees) to search across range [default = %(default)s]",
                      type=int, default=10)
    parser.add_argument("-w", "--wtranslations", dest="wtranslations",
                      help="Comma separated list of optimization weights of translations in "
                           "x, y, z for minctracc [default = %(default)s]",
                      type=str, default="0.2,0.2,0.2")
    parser.add_argument("--simplex", dest="simplex", type=float, default=1, 
                        help="Radius of minctracc simplex volume [default = %(default)s]")
    parser.set_defaults(use_multiple_seeds=True)
    parser.add_argument("--use-multiple-seeds", dest="use_multiple_seeds", action="store_true",
                        help="Find multiple possible starting points in the source and target for "
                             "the initial alignment in addition to using only the centre of gravity "
                             "(of the intensities) of the input files. [default = %(default)s]")
    parser.add_argument("--no-use-multiple-seeds", dest="use_multiple_seeds", action="store_false",
                        help="Opposite of --use-multiple-seeds")
    parser.set_defaults(max_number_seeds=3)
    parser.add_argument("--max-number-seeds", dest="max_number_seeds", type=int,
                        help="Specify the maximum number of seed-pair starting points "
                        "to use for the rotational part of the code. The seed "
                        "pairs are ordered based on the cross correlation gotten "
                        "from the alignment based on only the translation from the "
                        "seed point. [default = %(default)s]")
    parser.add_argument("source", help="", type=str, metavar="source.mnc")
    parser.add_argument("target", help="", type=str, metavar="target.mnc")
    parser.add_argument("output_xfm", help="", type=str, metavar="output.xfm")
    parser.add_argument("output_mnc", help="", type=str, metavar="output.mnc")

    options = parser.parse_args()

    if options.tmpdir:
        os.environ["TMPDIR"] = options.tmpdir
    elif "TMPDIR" not in os.environ:
        os.environ["TMPDIR"] = "/tmp/"

    print("TMP: %s" % os.environ["TMPDIR"])
    print("RANGE: %s INTERVAL: %s" % (options.range, options.interval))
    source     = options.source
    target     = options.target
    output_xfm = options.output_xfm
    output_mnc = options.output_mnc

    if options.resamplestepsize:
        source = downsample(options.source, options.resamplestepsize)
        target = downsample(options.target, options.resamplestepsize)
        # downsample the mask only if it is specified
        if options.mask:
            options.mask = downsample(options.mask, options.resamplestepsize)
    
    results = loop_rotations(stepsize=options.registrationstepsize, 
                             source=source,
                             target=target,
                             mask=options.mask, 
                             start=options.range,
                             interval=options.interval, 
                             wtranslations=options.wtranslations, 
                             simplex=options.simplex,
                             use_multiple_seeds=options.use_multiple_seeds,
                             max_number_seeds=options.max_number_seeds)
    
    print(results)
    subprocess.check_call(("cp %s %s" % (results[-1]["transform"], output_xfm)).split())
    subprocess.check_call(("cp %s %s" % (results[-1]["resampled"], output_mnc)).split())
    shutil.rmtree("%s/rot_%s" % (os.environ["TMPDIR"], os.getpid()))

