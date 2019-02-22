#!/usr/bin/env python3
#
# Convert between different mesh file formats
#
# author: jan.scholz@mouseimaging.ca
# original version: 2015-03-25
# based on: http://www.cfd-online.com/Forums/openfoam-meshing/86416-vtk-openfoam-something-else-can-reach-after-openfoam.html

from os import path
import sys
from optparse import OptionParser
import numpy as np
import vtk


class MyParser(OptionParser):
    """alow adding usage examples in epilog"""
    def format_epilog(self, formatter):
        return '\n' + self.epilog + '\n'


def view(stlfilename):
    reader = vtk.vtkSTLReader()
    reader.SetFileName(stlfilename)

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(reader.GetOutputPort())

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    # Create a rendering window and renderer
    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)

    # Create a renderwindowinteractor
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)

    # Assign actor to the renderer
    ren.AddActor(actor)

    # Enable user interface interactor
    iren.Initialize()
    renWin.Render()
    iren.Start()


def readColorFile(filename, verbose=False):
    ar = np.loadtxt(filename)
    if np.max(ar) > 255 or np.min(ar) < 0:
        raise ValueError('color values need to be in range [0..255]')
    if len(ar.shape) == 1:
        ar = np.vstack((ar,ar,ar)).T
    elif ar.shape[1] == 3:
        pass
    else:
        raise ValueError('use one or three column color array')
    return ar


def addColors(infilename, outfilename, colorfilename=None, colorstring=None,
                                                    binary=True, verbose=False):
    """add color array"""

    outformat = path.splitext(outfilename)[1].strip('.')
    if outformat!='ply':
        raise ValueError('colors are only supported for PLY format')

    informat = path.splitext(infilename)[1].strip('.')
    # set reader based on filename extension
    if informat=='stl':
        reader = vtk.vtkSTLReader()
    elif informat=='vtk':
        reader = vtk.vtkPolyDataReader()
    elif informat=='obj':
        reader = vtk.vtkMNIObjectReader()
    elif informat=='ply':
        reader = vtk.vtkPLYReader()
    elif informat=='vtp':
        reader = vtk.vtkXMLPolyDataReader()
    else:
        raise ValueError('cannot read input format: ' + informat)
    reader.SetFileName(infilename)
    reader.Update()

    #N = reader.GetOutput().GetNumberOfPolys()
    N = reader.GetOutput().GetNumberOfPoints()
    if verbose:
        print("read %i points (vertices)" % (N,))

    if colorfilename:
        colorar = readColorFile(colorfilename)
        if N != colorar.shape[0]:
            raise ValueError('number of rows in color file does not match' +
                             'number of points in mesh file')
    elif colorstring:
        color = [int(i) for i in colorstring.split()]
        colorar = np.ones((N,3)) * np.array(color)

    Colors = vtk.vtkUnsignedCharArray()
    Colors.SetNumberOfComponents(3)
    Colors.SetName("Colors")

    for i in range(0,N):
        Colors.InsertNextTuple3(*colorar[i,:])

    polydata = vtk.vtkPolyData()
    polydata = reader.GetOutput()

    polydata.GetPointData().SetScalars(Colors)
    polydata.Modified()

    writer = vtk.vtkPLYWriter()
    writer.SetArrayName("Colors")
    writer.SetInputData(polydata)
    writer.SetFileName(outfilename)
    if binary:
        if verbose: print('setting output to binary')
        writer.SetFileTypeToBinary()
    else:
        if verbose: print('setting output to ascii')
        writer.SetFileTypeToASCII()
    err = writer.Write()


def readMeshFile(filename, clean=True, verbose=False, recompute_normals=True):
    """Read mesh file.
    The input format is determined by file name extension.
    Polygons get split into triangles to support various restrictive output
    formats.
    If clean, degenerate data gets removed."""

    informat = path.splitext(options.infilename)[1].strip('.')
    # set reader based on filename extension
    if informat=='stl':
        reader = vtk.vtkSTLReader()
    elif informat=='vtk':
        reader = vtk.vtkPolyDataReader()
    elif informat=='obj':
        reader = vtk.vtkMNIObjectReader()
    elif informat=='ply':
        reader = vtk.vtkPLYReader()
    elif informat=='vtp':
        reader = vtk.vtkXMLPolyDataReader()
    #elif informat=='tag':
    #    reader = vtk.vtkMNITagPointReader()
    else:
        raise ValueError('cannot read input format: ' + informat)
    reader.SetFileName(filename)
    reader.Update()

    if verbose:
        print("read %i polygons from file %s" % \
                               (reader.GetOutput().GetNumberOfPolys(), filename))

    # merge duplicate points, and/or remove unused points and/or remove degenerate cells
    if clean:
        polydata = vtk.vtkCleanPolyData()
        polydata.SetInputConnection(reader.GetOutputPort())
        poly_data_algo = polydata
        if verbose:
            print("cleaned poly data")
    else:
        poly_data_algo = reader

    # convert input polygons and strips to triangles
    triangles = vtk.vtkTriangleFilter()
    triangles.SetInputConnection(poly_data_algo.GetOutputPort())

    if recompute_normals:
        normals = vtk.vtkPolyDataNormals()
        normals.SetInputConnection(triangles.GetOutputPort())

        normals.SplittingOff()
        normals.ComputePointNormalsOn()
        normals.AutoOrientNormalsOn()
        normals.ConsistencyOn()
        normals.NonManifoldTraversalOn()
        if verbose:
            print("recomputed normals")
            print("finished reading", filename)
        return normals
    else:
        if verbose:
            print("finished reading", filename)
        return triangles


def writeMeshFile(triangles, filename, binary=True, verbose=False):
    """Write mesh file.
    The output format is determined by file name extension. Files can be written
    in binary (default) and ASCII format."""

    outformat = path.splitext(options.outfilename)[1].strip('.')
    # set writer based on filename extension
    if outformat=='stl':
        write = vtk.vtkSTLWriter()
    elif outformat=='vtk':
        write = vtk.vtkPolyDataWriter()
    elif outformat=='obj':
        write = vtk.vtkMNIObjectWriter()
    elif outformat=='ply':
        write = vtk.vtkPLYWriter()
    elif outformat=='vtp':
        write = vtk.vtkXMLPolyDataWriter()
    elif outformat=='tag':
        write = vtk.vtkMNITagPointWriter()
    else:
        raise ValueError('cannot write output format' + outformat)
    write.SetInputConnection(triangles.GetOutputPort())

    if outformat!='tag':
        if binary:
            if verbose: print('setting ouptut to binary')
            write.SetFileTypeToBinary()
        else:
            if verbose: print('setting ouptut to ascii')
            write.SetFileTypeToASCII()

    write.SetFileName(filename)
    err = write.Write()
    if err != 1:
        raise IOError('failed to write')

    if verbose:
        print("wrote", filename)
    pass



if __name__ == "__main__":
    usage = """usage: %prog [-h/--help] -i INFILE -o OUTFILE"""

    description = """Convert between mesh file formats.
    Currently supports reading and writing of STL, VTK, OBJ (BIC object)"""
    epilog = "Example:\n  " + \
        path.basename(__file__) + " -v --ascii -i foo.vtk -o bar.stl"

    parser = MyParser(usage=usage, description=description, epilog=epilog)

    parser.add_option("-i", "--input", dest="infilename",
                      help="no help",
                      type='string', default="")
    parser.add_option("-o", "--output", dest="outfilename",
                      help="no help",
                      type='string', default="")
    parser.add_option("-a", "--ascii", dest="binary",
                      help="save in ascii format",
                      action="store_false", default=True)
    parser.add_option("--noclean", dest="clean",
                      help="remove degenerate data",
                      action="store_false", default=True)
    parser.add_option("--no-recompute-normals", dest="recompute_normals",
                      help="recompute surface normals",
                      action="store_false", default=True)
    parser.add_option("--color", dest="color",
                      help="add color, three values [0..255]",
                      type='string', default=None)
    parser.add_option("--colorfile", dest="colorfile",
                      help="add color, three column file [0..255]",
                      type='string', default=None)
    parser.add_option("--view", dest="view",
                      help="view stl file",
                      action="store_true", default=False)
    parser.add_option("-v", "--verbose", dest="verbose",
                      help="more verbose output",
                      action="store_true", default=False)

    (options, args) = parser.parse_args()

    if options.infilename is '':
        parser.error('INFILE not specified (-i)')

    if options.outfilename is '':
        parser.error('OUTFILE not specified (-o)')

    if not path.exists(options.infilename):
        parser.error('could not find INFILE')

    outpath = path.dirname(options.outfilename)
    if outpath and not path.exists(outpath):
        parser.error('output directory does not exist: ' + outpath)

    if options.view:
        view(options.infilename)
        sys.exit(0)

    if options.color or options.colorfile:
        addColors(options.infilename, options.outfilename,
                  colorstring=options.color,
                  colorfilename=options.colorfile,
                  binary=options.binary, verbose=options.verbose)
        sys.exit(0)

    triangles = readMeshFile(options.infilename, clean=options.clean, verbose=options.verbose, recompute_normals=options.recompute_normals)

    writeMeshFile(triangles, options.outfilename, binary=options.binary,
                  verbose=options.verbose)
