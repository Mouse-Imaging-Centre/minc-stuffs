#!/usr/bin/env python3

import argparse, subprocess, os
from typing import Tuple

def run_subprocess(cmds):
    cmdstr = " ".join(cmds)
    print("\n")
    print(cmdstr)
    p = subprocess.Popen(cmdstr, shell=True)
    if p.wait() != 0:
        raise Exception("Something went wrong!")
    return p

def explode(filename: str) -> Tuple[str, str, str]:
    base, ext = os.path.splitext(filename)
    directory, name = os.path.split(base)
    return (directory, name, ext)

if __name__ == "__main__":

    description = 'Compute the determinant of a transform.'

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('--clobber',
                        action='store_true',
                        dest='clobber',
                        help='clobber all temporary and output files'
                        )
    parser.add_argument('--inverse',
                        action='store_true',
                        dest='inverse',
                        help='compute the determinant of the inverse instead'
                        )
    parser.add_argument('--non-linear-only',
                        action='store_true',
                        default=False,
                        dest='non_linear_only',
                        help='compute the determinant of the non-linear portion only. please provide a mask using --mask'
                        )
    parser.add_argument('--mask',
                        nargs="?",
                        dest='mask',
                        help='mask for computing the non-linear portion of the input transform'
                        )
    parser.add_argument('--smooth',
                        dest='smooth',
                        help='smooth the displacement field before calculating the determinant'
                        )
    parser.add_argument('--log',
                        action='store_true',
                        default=False,
                        dest='log',
                        help='compute the logarithm of the transformation'
                        )
    parser.add_argument("--temp-folder",
                        dest="temp_folder",
                        default="/tmp")

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("--like",
                               dest="input_like",
                               help='a like file',
                               required=True)
    requiredNamed.add_argument("--transform",
                               dest="input_transform",
                               help='a transform for which the determinant should be calculated',
                               required=True)
    requiredNamed.add_argument("--determinant",
                               dest='output_determinant',
                               help='output',
                               required=True)

    args = parser.parse_args()

    tempdir = args.temp_folder
    if not os.path.isdir(tempdir):
        os.mkdir(tempdir)

    input_dir, input_name, input_ext = explode(args.input_transform)
    output_dir, output_name, output_ext = explode(args.output_determinant)

# These 4 cases are the product of 2 optional preprocessing steps (i know its ugly)
    if args.inverse and args.non_linear_only:
        inverse = os.path.join(tempdir, input_name + "_inverse.xfm")
        p = run_subprocess(["xfminvert",
                            "-clobber" if args.clobber else "",
                            args.input_transform,
                            inverse])
        linear_part = os.path.join(tempdir, input_name + "_linear_part.xfm")
        p = run_subprocess(["lin_from_nlin", "-lsq12",
                            "-clobber" if args.clobber else "",
                            "-mask", args.mask,
                            args.input_like,
                            args.input_transform,
                            linear_part])
        transform = os.path.join(tempdir, input_name + "_nlin_part.xfm")
        p = run_subprocess(["xfmconcat",
                            "-clobber" if args.clobber else "",
                            inverse,
                            linear_part,
                            transform])
    elif args.non_linear_only:
        linear_part = os.path.join(tempdir, input_name + "_linear_part.xfm")
        p = run_subprocess(["lin_from_nlin", "-lsq12",
                            "-clobber" if args.clobber else "",
                            "-mask", args.mask,
                            args.input_like,
                            args.input_transform,
                            linear_part])

        linear_part_inv = os.path.join(tempdir, input_name + "_linear_part_inv.xfm")
        p = run_subprocess(["xfminvert",
                            "-clobber" if args.clobber else "",
                            linear_part,
                            linear_part_inv])
        transform = os.path.join(tempdir, input_name + "_nlin_part.xfm")
        p = run_subprocess(["xfmconcat",
                            "-clobber" if args.clobber else "",
                            linear_part_inv,
                            args.input_transform,
                            transform])
    elif args.inverse:
        transform = os.path.join(tempdir, input_name + "_inverse.xfm")
        p = run_subprocess(["xfminvert",
                            "-clobber" if args.clobber else "",
                            args.input_transform,
                            transform])
    else:
        transform = args.input_transform


# Start the linear pipeline
    displacement = os.path.join(tempdir, output_name + "_displacement.mnc")
    p = run_subprocess(["minc_displacement",
                        "-clobber" if args.clobber else "",
                        args.input_like,
                        transform,
                        displacement])

    if args.smooth:
        smooth_displacement = os.path.join(tempdir, output_name + "_smooth_displacement.mnc")
        p = run_subprocess(["smooth_vector",
                            "--clobber" if args.clobber else "",
                            "--filter", "--fwhm", args.smooth,
                            displacement,
                            smooth_displacement])

    temp_determinant = os.path.join(tempdir, output_name + "_temp_determinant.mnc")
    p = run_subprocess(["mincblob",
                        "-determinant",
                        "-clobber" if args.clobber else "",
                        smooth_displacement if args.smooth else displacement,
                        temp_determinant])
    if not args.log:
        p = run_subprocess(["mincmath -2 -const 1 -add", #mincblob output needs this
                                 "-clobber" if args.clobber else "",
                            temp_determinant,
                            args.output_determinant
                            ])
    else:
        nolog_determinant = os.path.join(tempdir, output_name + "_nolog_determinant.mnc")
        p = run_subprocess(["mincmath -2 -const 1 -add",  # mincblob output needs this
                            "-clobber" if args.clobber else "",
                            temp_determinant,
                            nolog_determinant
                            ])
        p = run_subprocess(["mincmath -2 -log",
                            "-clobber" if args.clobber else "",
                            nolog_determinant,
                            args.output_determinant
                            ])