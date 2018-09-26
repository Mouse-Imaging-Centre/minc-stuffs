#!/usr/bin/env python3

import argparse, subprocess, os

def run_subprocess(cmds):
    cmdstr = " ".join(cmds)
    print("\n")
    print(cmdstr)
    p = subprocess.Popen(cmdstr,
                         shell=True)
    p.wait()
    return p
    # (out,err)=p.communicate()
    # out=out.decode()
    # if p.wait() != 0:
    #     print(err.decode())

if __name__ == "__main__":

    usage = "%(prog)s [options] input_transform.xfm input_like.mnc output_determinant.mnc"
    description = 'Compute the determinant of a transform.'

    parser = argparse.ArgumentParser(usage=usage, description=description)

    parser.add_argument('input_like',
                        help='a like file',
                        )
    parser.add_argument('input_transform',
                        type=str,
                        help='a transform for which the determinant should be calculated',
                        )
    parser.add_argument('output_determinant',
                        help='output',
                        )
    parser.add_argument('--clobber',
                        action='store_true',
                        dest='clobber',
                        help='clobber all temporary and output files'
                        )
    parser.add_argument('--non-linear-only',
                        action='store_true',
                        default=False,
                        dest='non_linear_only',
                        help='compute the determinant of the non-linear portion only'
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

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('--temp-folder',
                               dest ="temp_folder",
                               help='temp folder',
                               required=True)

    args = parser.parse_args()

    tempdir = args.temp_folder
    if not os.path.isdir(tempdir):
        os.mkdir(tempdir)

    no_ext, ext = os.path.splitext(args.output_determinant)
    output_dir, output_name = os.path.split(no_ext)

###########START
    #if non_linear_only:


    displacement = os.path.join(tempdir, output_name + "_displacement.mnc")
    p = run_subprocess(["minc_displacement",
                        "-clobber" if args.clobber else "",
                        args.input_like,
                        args.input_transform,
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