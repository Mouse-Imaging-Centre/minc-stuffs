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
    line = lines.split('\n')[-2]
    cog = array(line.strip().split(" ")).astype("float32")
    return cog

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

def minctracc(source, target, mask, stepsize, wtranslations):
    wtrans_decomp = array(wtranslations.split(',')).astype("float")
    tmp_transform = get_tempfile('.xfm')
    if mask is not None:
        cmd = ("minctracc -identity -lsq6 -xcorr -step %s %s %s %s %s %s -source_mask %s -model_mask %s -w_translations %s %s %s"
                % (stepsize, stepsize, stepsize, source, target, tmp_transform, mask, mask,
                   wtrans_decomp[0], wtrans_decomp[1], wtrans_decomp[2]))
        print(cmd)
        subprocess.check_call(cmd.split())
    else:
        cmd = ("minctracc -identity -lsq6 -xcorr -step %s %s %s %s %s %s -w_translations %s %s %s"
                 % (stepsize, stepsize, stepsize, source, target, tmp_transform,
                    wtrans_decomp[0], wtrans_decomp[1], wtrans_decomp[2]))
        print(cmd)
        subprocess.check_call(cmd.split())
    return tmp_transform

def concat_transforms(t1, t2):
    tmp_transform = get_tempfile('.xfm')
    subprocess.check_call(("xfmconcat %s %s %s" % (t1, t2, tmp_transform)).split())
    return tmp_transform
              
def loop_rotations(stepsize, source, target, mask, start=50, interval=10, wtranslations="0.2,0.2,0.2"):

    # get the centre of gravity for both volumes
    cog1 = get_centre_of_gravity(source)
    cog2 = get_centre_of_gravity(target)
    cogdiff = cog2 - cog1
    print("\n\nCOG diff: %s\n\n" % cogdiff.tolist())

    # load the target and mask volumes
    targetvol = volumeFromFile(target)
    maskvol = volumeFromFile(mask) if mask is not None else None
    
    results = []
    best_xcorr = 0
    for x in range(-start, start+1, interval):
        for y in range(-start, start+1, interval):
            for z in range(-start, start+1, interval):
                # we need to include the centre of the volume as rotation centre = cog1
                init_transform = create_transform(cogdiff, x, y, z, cog1)
                #init_resampled = resample_volume(source,target, init_transform)
                init_resampled = resample_volume(source ,target, init_transform)
                transform = minctracc(init_resampled, target, mask, stepsize, wtranslations)
                resampled = resample_volume(init_resampled, target, transform)
                xcorr = compute_xcorr(resampled, targetvol, maskvol)
                if isnan(xcorr):
                    xcorr = 0
                conc_transform = concat_transforms(init_transform, transform)
                results.append({'xcorr': xcorr, 'transform': conc_transform, \
                                'resampled': resampled, 'x': x, \
                                'y': y, 'z': z})
                if xcorr > best_xcorr:
                    best_xcorr = xcorr
                else:
                    os.remove(resampled)
                os.remove(init_resampled)
                print("FINISHED: %s %s %s :: %s" % (x,y,z, xcorr))
    targetvol.closeVolume()
    if mask is not None:
        maskvol.closeVolume()
    sort_results(results)
    return results

def dict_extract(adict, akey):
    return adict[akey]

def sort_results(results):
    sort_key_func = functools.partial(dict_extract, akey="xcorr")
    results.sort(key=sort_key_func)

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
                      help="resampled volumes to this stepsize",
                      type=float, default=0.2)
    parser.add_argument("-g", "--stepsizeregistration", dest="registrationstepsize",
                      help="use this stepsize in the minctracc registration",
                      type=float, default=0.6)
    parser.add_argument("-t", "--tempdir", dest="tmpdir",
                      help="temporary directory to use",
                      type=str)
    parser.add_argument("-r", "--range", dest="range",
                      help="range of rotations to search across",
                      type=int, default=50)
    parser.add_argument("-i", "--interval", dest="interval",
                      help="interval (in degrees) to search across range",
                      type=int, default=10)
    parser.add_argument("-w", "--wtranslations", dest="wtranslations",
                      help="Comma separated list of optimization weights of translations in x, y, z for minctracc",
                      type=str, default="0.2,0.2,0.2")
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
    if options.mask:
        results = loop_rotations(options.registrationstepsize, source, target, options.mask, options.range,
                                 options.interval, options.wtranslations)
    else:
        results = loop_rotations(options.registrationstepsize, source, target, None, options.range,
                                 options.interval, options.wtranslations)
    
    print(results)
    subprocess.check_call(("cp %s %s" % (results[-1]["transform"], output_xfm)).split())
    subprocess.check_call(("cp %s %s" % (results[-1]["resampled"], output_mnc)).split())
    shutil.rmtree("%s/rot_%s" % (os.environ["TMPDIR"], os.getpid()))

