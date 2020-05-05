"""
**LineStacker Module:**\n
Main Stacker module. Contains all basic functions.
"""
# -*- coding: utf-8; -*-
# stacker, Python module for stacking of interferometric data.
# Copyright (C) 2014  Lukas Lindroos
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor,
# Boston, MA 02110-1301, USA.

"""
    Library to stack interferometric images.
"""

import math
import os
from ctypes import cdll
import re
import glob
import numpy as np

__author__ = 'Lukas Lindroos'
__copyright__ = 'Copyright 2014'
__license__ = 'GPL'
__version__ = '1.0.3'
__maintainer__ = 'Lukas Lindroos'
__email__ = 'lindroos@chalmers.se'

PB_CONST = 0
PB_MS = 1
PB_FITS = 2

FILETYPENAME = {}
FILE_TYPE_NONE = 0
FILETYPENAME[FILE_TYPE_NONE] = 'none'
FILE_TYPE_MS = 1
FILETYPENAME[FILE_TYPE_MS] = 'ms'
FILE_TYPE_FITS = 2
FILETYPENAME[FILE_TYPE_FITS] = 'fits'

MS_DATACOLUMN_DATA = 1
MS_MODELCOLUMN_DATA = 2


clib_path = os.path.join(os.path.abspath(__path__[0]),
                         'stacker_clib')


def _casa_svnversion(casapath):
    import re
    import os.path
    raise Exception
    try:
        casapyinfofile = open(os.path.join(casapath, 'casapyinfo'))
        svnversion = None
        for line in casapyinfofile:
            match = re.match('SVNVERSION="([0-9]*)"', line)
            if match:
                svnversion = match.groups()[0]
    except IOError:
        casaconfigfile = open(os.path.join(casapath, 'bin', 'casa-config'))
        for line in casaconfigfile:
            match = re.match('.*\$revision="([0-9]*)";.*', line)
            if match:
                svnversion = match.groups()[0]

    if svnversion is None:
        raise IOError('Can not find casa version.')

    return svnversion


def _libs_matching_svn_version(svnversion):
    """Return the list of stacker versions for the specified svn version"""
    return sorted(glob.glob(os.path.join(clib_path, 'libstacker-r{0}*.so'.format(svnversion))), reverse=True)


def _load_stacker_casa_svn_version(svnversion):
    """Sequentially attempt to load the options for the supplied svn version until successful"""
    libpath_options = _libs_matching_svn_version(svnversion)
    print('Loading stacking library for casapy svn revision {0}'.format(svnversion))
    for libpath in libpath_options:
        try:
            libstacker = cdll.LoadLibrary(libpath)
            print("{0} loaded".format(libpath))
            return libstacker
        except OSError:
            continue
    else:
        print("Could not load library for svn version {0}".format(svnversion))

def _load_stacker_lib():
    """Attempt to load a suitable stacker lib"""
    try:
        # If we are inside CASA - use libs linked against NRAO's CASA releases
        from taskinit import casa
        svnversion = _casa_svnversion(casa['dirs']['root'])
        lib_list = _libs_matching_svn_version(svnversion)
        if len(lib_list) > 0:
            libstacker = _load_stacker_casa_svn_version(svnversion)
        else:
            print('warning, no precompiled library compatible with your version of casa exists.')
            print('It is recommended to recompile stacker for your casa version.')
            stacker_clib_ls = os.listdir(clib_path)
            stackerlibs = [((re.match('libstacker-r([0-9]*).*.so', f).group(1)), f)
                           for f in stacker_clib_ls
                           if re.match('libstacker-r[0-9]*.so', f)]
            print(stackerlibs)
            vdiff = [(abs(int(svnversion)-int(v)), lib, v)
                     for (v, lib) in stackerlibs]
            vdiff.sort()
            print(vdiff)
            print('Trying to use svn revision {0}.'.format(vdiff[0][2]))
            libpath = os.path.join(clib_path, vdiff[0][1])
            try:
                libstacker = cdll.LoadLibrary(libpath)
            except OSError:#, e:
                print(e)
                print("Loading libstacker failed. You may need to build stacker for your version of CASA.")
                return
    except ImportError:
        # We are in a pure python session.
        # Not that there is anything wrong with that.
        libpath = os.path.join(clib_path, 'libstacker.so')
        libstacker = cdll.LoadLibrary(libpath)
    return libstacker


#libstacker = _load_stacker_lib()


class CoordList(list):
    """
        Extended list to contain list of coordinates.

        Requires an image list in case of pixel coordinates to work properly.

        Parameters
        ---------
        imagenames
            A list of image names, requiered for pixle coordinates to work properly.
        coord_type
            **'physical'** or **'pixel'**. **'physical'** coordinates should be converted to pixels using **LineStacker.getPixelCoords** before stacking.
        unit
            can be **'rad'**, **'deg'** or **'pix'**.

    """
    def __init__(   self,
                    imagenames=[],
                    coord_type='physical',
                    unit='rad'):

        super(CoordList, self).__init__()

        if isinstance(imagenames, str):
            imagenames = [imagenames]

        self.coords = []
        self.imagenames = imagenames
        self.coord_type = coord_type
        self.unit = unit

    def __getitem__(self, i):
        return self.coords[i]

    def __setitem__(self, i, x):
        self.coords[i] = x

    def append(self, x):
        self.coords.append(x)

    def __len__(self):
        return len(self.coords)

    def __iter__(self):
        for x in self.coords:
            yield x

    def __getslice__(self, i, j):
        new_a = CoordList(self.imagenames, self.coord_type, self.unit)
        new_a.coords = self.coords.__getslice__(i, j)
        return new_a

    def __repr__(self):
        ret = []
        for x in self.coords:
            ret.append('(' + x.__str__() + ')')
        return '\n'.join(ret)

    def __str__(self):
        ret = []
        for x in self.coords:
            ret.append('(' + x.__str__() + ')')
        return '\n'.join(ret)
#         return '{0}, {1}'.format(self.x, self.y)


class Coord:
    """
        Describes a stacking position on an image.

        Class used internally to represent coordinates. May describe a
        physical coordinate or a pixel coordinate.

        Init: Creates a coordinate. A pixel coordinate should always specify
        to which image it belongs. Physical coordinates should be in
        J2000 radians.

        Parameters
        ---------
        x
            x coordinate of stacking target
        y
            y coordinate of stacking target
        z
            redshift of stacking target
        obsSpecArg
            argument of spectral bin on which to center the stack, will be computed from **z** and **fEm** if not specified
        weight
            Weight of the. source in case of mean stacking.
            Default is 1.
        image
            Index of the image associated to the source.
            Automatically set with **LineStacker.readCoords**.

    """

    def __init__(self, x, y, z=0, obsSpecArg=0, weight=1, image=0):
        self.x = x
        self.y = y
        self.z = z
        self.weight = weight
        self.image = image
        self.obsSpecArg=obsSpecArg

    def __str__(self):
        return '{0}, {1}'.format(self.x, self.y)

    def setWeightToX(self, X, chanWidth=0):
        if chanWidth==0:
            self.weight=X
        else:
            self.weight=np.array([X for i in range(chanWidth)])

    def setZeroWeightLeft(self,numberOfZeroLeftF, freqlen):
        try:
            len(self.weight)
            if numberOfZeroLeftF!=0:
                self.weight=np.concatenate((np.zeros(abs(numberOfZeroLeftF)), self.weight ))
        except TypeError: #if weights are not list
            self.weight=np.ones(freqlen)*self.weight
            if numberOfZeroLeftF!=0:
                #self.weight[0:abs(numberOfZeroLeftF)]=0
                self.weight=np.concatenate((np.zeros(abs(numberOfZeroLeftF)), self.weight ))

    def setZeroWeightRight(self,numberOfZeroRightF, freqlen):
        try:
            len(self.weight)
            if numberOfZeroRightF!=0:
                self.weight=np.concatenate((self.weight, np.zeros(numberOfZeroRightF)  ))

        except TypeError: #if weights are not list
            self.weight=np.ones(freqlen)*self.weight
            if numberOfZeroRightF!=0:
            #self.weight[-numberOfZeroRightF:]=0
                self.weight=np.concatenate((self.weight, np.zeros(numberOfZeroRightF)  ))

def readCoordsGUI(unit='deg', lineON=True, **kwargs):
    import Tkinter, Tkconstants, tkFileDialog
    filez=  tkFileDialog.askopenfilenames(initialdir = ".",title = "Select coord files",filetypes = (("txt files","*.txt"),("all files","*.*")))
    filez=list(filez)

    coords=readCoords(filez, unit=unit, lineON=lineON, **kwargs)
    return coords

def readCoordsNamesGUI():
    """
        Open GUI to select coordinates files.

        Returns path of selected files.
    """
    import Tkinter, Tkconstants, tkFileDialog
    filez=  tkFileDialog.askopenfilenames(initialdir = ".",title = "Select coord files",filetypes = (("txt files","*.txt"),("all files","*.*")))
    filez=list(filez)
    return filez


#def selectImagesGUI():
#    import Tkinter, Tkconstants, tkFileDialog
#    filez=  tkFileDialog.askopenfilenames(initialdir = ".",title = "Select image files",filetypes = (("casa im cubes","*.im"),("casa image cubes","*.image"),("all files","*.*")))
#    return list(filez)

def readCoords(coordfiles, unit='deg', lineON=True, delimiter='default'):
    if delimiter=='default':
        delimiter=[',',' ','\t']
        delimiterNames=['coma', 'space', 'tab']
    else:
        delimiter=[delimiter]
    for (j,delimiter) in enumerate(delimiter):

        try:
            coords=_readCoords(coordfiles, unit=unit, lineON=True, delimiter=delimiter)
            return coords
        except RuntimeError:
            print "\033[1;33m"+"\nUsing delimiter '"+str(delimiterNames[j])+"' lead to an error, trying other delimiter (if any)"+"\x1b[0m"+'\n'
            pass
    raise Exception("RuntimeError, check data format. If set to default delimiters can be comma, space, and tab. Otherwise feed own delimiter to this function.")

def _readCoords(coordfiles, unit='deg', lineON=True, delimiter=','):
    """
        Reads coordinate files from disk and produces a list./!\\\\ To each image should be associated one single coord file.

        Parameters
        ---------
        coordfiles:
            List of path to coordinate files. Files in csv format. x and y should be in J2000. If using to stack line, redshift should be in the third column. A weight may be added in the last column to weight positions for mean stacking. If set to **None** stacking position is defined as the center of each cube.
        unit:
            Unit of input coordinates. Allows two values, **'deg'** and **'rad'**.
        lineON
            Either **True** or **False**, should be set to **True** if doing line stacking. Default is **True**.
    """
    import csv

    coords = CoordList()
    if coordfiles==None:
        print 'coordfile is None, stacking on the center'
        return None
    for (i, coordi) in enumerate(coordfiles):
        coordreader = csv.reader(open(coordi, 'rb'), delimiter=delimiter)
        for row in coordreader:
            try:
                x = float(row[0])
                y = float(row[1])
            except Exception:
                from taskinit import qa
                x=qa.convert(row[0], 'rad')['value']
                y=qa.convert(row[1], 'rad')['value']
                unit='rad'
            if unit == 'deg':
                x = math.pi/180.*x
                y = math.pi/180.*y

            if x > 2*math.pi:
                x -= 2*math.pi
            if y > math.pi:
                y -= 2*math.pi

            if not lineON:
                if len(row)==3:
                    weight = float(row[2])
                    coords.append(Coord(x, y, weight=weight, image=i))
                elif len(row)==2:
                    weight = 1.
                    coords.append(Coord(x, y, weight=weight, image=i))
                elif len(row)==4:
                    lineON=True
                    print '/!\\ care : Line mode is set to False but 4 coordinate rows are found, line mode is now set to True'
                else:
                    raise Exception('too many rows in coordfile, structure should be x, y, z, w')

            elif lineON:
                if len(row) == 4:
                    weight = float(row[3])
                elif len(row)==3:
                    weight = 1.
                else:
                    raise Exception('too many rows in coordfile, structure should be x, y, z, w')
                z = float(row[2])
                coords.append(Coord(x, y, z=z, weight=weight, image=i))

    return coords


def writeCoords(coordpath, coords, unit='deg'):
    """
        Write coordinates to a file

        Parameters
        ---------
        coordpath
            absolute path to coordinate file to be created
        coords
            list of coordinates (stacker.Coord instances) to write to file
    """
    import csv

    with open(coordpath, 'wb') as coordfile:
        coordwriter = csv.writer(coordfile, delimiter=',')
        for coord in coords:
            x = coord.x
            y = coord.y
            z=coord.z
            weight = coord.weight

            if unit == 'deg':
                x = x*180./math.pi
                y = y*180./math.pi

            coordwriter.writerow([str(x), str(y), str(z), str(weight)])


def _checkfile(filename, datacolumn):
    import re
    # Currently this supports only ms files
    # As such there is no reason to check filetype.
    # If it cannot be opened as ms it will not be supported.
#     if re.match('^.*[mM][sS]/*$', filename) is not None:
    try:
        from taskinit import ms
        ms.open(filename)
        ms.done()
    except ImportError:
        # This probably means that it was run from a pure python session.
        # We will relegate any checks that it is a valid ms file to the
        # stacker.
        if not os.access(filename, os.F_OK):
            raise RuntimeError('Could not find data file "{}".'.format(
                filename))

    filename = filename
    filetype = FILE_TYPE_MS
    fileoptions = 0
    if datacolumn == 'data':
        fileoptions = MS_DATACOLUMN_DATA
    elif datacolumn == 'model' or datacolumn == 'model_data':
        fileoptions = MS_MODELCOLUMN_DATA
#     elif re.match('^.*[fF][iI][tT][sS]$', filename) is not None:
#         raise NotImplementedError('FITS format is currently not supported.')
    return filetype, filename, fileoptions

def coordsTocl(name, flux, coords):
    from taskinit import cl, qa

    flux = qa.quantity(flux)
    cl.done()
    cl.rename(name)
    for coord in coords:
        clpars = {}
        clpars['flux'] = -flux['value']
        clpars['fluxunit'] = flux['unit']
        clpars['dir'] = ['J2000', str(coord.x)+'rad', str(coord.y)+'rad']
        clpars['shape'] = 'point'
        cl.addcomponent(**clpars)
    cl.done()


def randomCoords(imagenames, ncoords=10):
    """
        Randomize a set of coordinates anywhere on any images.

        Parameters
        ---------
        imagenames
            A list of images paths.
        ncoords
            Number of random coordinates. Default is 10.
    """
    import random
    from taskinit import ia, qa

    xmin, xmax = [], []
    ymin, ymax = [], []
    minObsSpecArg, maxObsSpecArg=[],[]
    for image in imagenames:
        ia.open(image)
        trc = ia.boundingbox()['trcf'].split(', ')
        blc = ia.boundingbox()['blcf'].split(', ')
        xmin.append(qa.convert(qa.quantity(trc[0]), 'rad')['value'])
        xmax.append(qa.convert(qa.quantity(blc[0]), 'rad')['value'])
        ymin.append(qa.convert(qa.quantity(blc[1]), 'rad')['value'])
        ymax.append(qa.convert(qa.quantity(trc[1]), 'rad')['value'])
        minObsSpecArg.append(0)
        maxObsSpecArg.append(ia.boundingbox()['regionShape'][-1]-1)
        ia.done()

    randomcoords = CoordList(imagenames)
    for i in range(ncoords):
        imageid = random.randint(0, len(imagenames)-1)
        x = random.uniform(xmin[imageid], xmax[imageid])
        y = random.uniform(ymin[imageid], ymax[imageid])
        obsSpecArg=random.randint(minObsSpecArg[imageid], maxObsSpecArg[imageid])
        c = Coord(x, y, obsSpecArg=obsSpecArg)
        randomcoords.append(c)

    return randomcoords


def randomizeCoords(coords, beam,maxBeamRange=5):
    """
        Randomize a new set of coordinates at a distance
        [beam, maxBeamRange*beam] of the original coordinates

        Parameters
        ---------
        coords
            list of original coordinates (stacker.Coord instances)
        beam
            beam size is radians,
            new random coordinates will be at a minimum distance beam
            from the original coordinates.
        maxBeamRange
            maximum distance from original coordinates
            at which new coordinates can be located, in units of beams.
            Default is 5.
    """

    import random
    import math

    xmin, xmax = [], []
    ymin, ymax = [], []
    minObsSpecArg, maxObsSpecArg=[],[]

    for image in coords.imagenames:
        ia.open(image)
        trc = ia.boundingbox()['trcf'].split(', ')
        blc = ia.boundingbox()['blcf'].split(', ')
        xmin.append(qa.convert(qa.quantity(trc[0]), 'rad')['value'])
        xmax.append(qa.convert(qa.quantity(blc[0]), 'rad')['value'])
        ymin.append(qa.convert(qa.quantity(blc[1]), 'rad')['value'])
        ymax.append(qa.convert(qa.quantity(trc[1]), 'rad')['value'])
        minObsSpecArg.append(0)
        maxObsSpecArg.append(ia.boundingbox()['regionShape'][-1]-1)
        ia.done()

    randomcoords = CoordList(coords.imagenames, coords.coord_type,
                             unit=coords.unit)

    for coord in coords:
        dr = random.uniform(beam, maxBeamRange*beam)
        dphi = random.uniform(0, 2*math.pi)
        x = coord.x + dr*math.cos(dphi)
        y = coord.y + dr*math.sin(dphi)
        tempCounter=0
        while (     x<xmin[coord.image] or x>xmax[coord.image] or
                    y<ymin[coord.image] or y>ymax[coord.image]):
            dr = random.uniform(beam, 5*beam)
            dphi = random.uniform(0, 2*math.pi)
            x = coord.x + dr*math.cos(dphi)
            y = coord.y + dr*math.sin(dphi)
            tempCounter+=1
            if tempCounter>100:
                raise Exception('Cannot find random position, change search region')
        obsSpecArg=random.randint(minObsSpecArg[coord.image], maxObsSpecArg[coord.image])
        randomcoords.append(Coord(x, y, obsSpecArg=obsSpecArg, weigth=coord.weight, image=coord.image))

    return randomcoords



def make_pbfile(vis, pbfile):
    from taskinit import im, ms, ia, qa, tb
    import numpy as np
    from scipy.constants import c

    ms.open(vis)
    fields = ms.range('field_id')['field_id']
    ms.done()
    im.open(vis)
    im.selectvis(field=fields[0])
    ms.open(vis)
    freq = np.mean(ms.range('chan_freq')['chan_freq'])
    phase_dir = ms.range('phase_dir')['phase_dir']['direction']
    ms.done()

    phase_dir = phase_dir[0][0], phase_dir[1][0]
    phase_dir = [qa.formxxx(str(phase_dir[0])+'rad', format='hms'),
                 qa.formxxx(str(phase_dir[1])+'rad', format='dms')]
    phase_dir = 'J2000 '+' '.join(phase_dir)

    tb.open(vis+'/ANTENNA/')
    dishdia = np.min(tb.getcol('DISH_DIAMETER'))
    tb.done()

    # pb of 512 pix cover pb down to 0.001
    # ensure largest pixel to pixel var to .01
    minpb = 0.001
    nx = 512
    cellconv = (nx*np.sqrt(np.log(2)/np.log(1/minpb)))**-1

    beam = c/freq/dishdia
    cell = {}
    cell['value'] = beam*cellconv
    cell['unit'] = 'rad'

#     nx = int(3*3e8/freq/dishdia*1.22*180/
#              math.pi*3600/qa.convert(advise['cell'],
#              'arcsec')['value'])
    # Chosen as to be 3 times fwhm of primary beam,
    # should include up to approximately .01 of peak flux

    im.defineimage(nx=nx, ny=nx, cellx=cell, celly=cell, phasecenter=phase_dir)
    im.setvp(dovp=True)
    im.makeimage(type='pb', image=pbfile)
    im.done()
    ia.open(pbfile)
    cs = ia.coordsys()
    cs.setreferencevalue(type='direction', value=[0., 0.])
    ia.setcoordsys(cs.torecord())
    ia.maskhandler('delete', 'mask0')
    ia.done()


def getPixelCoords(coords, imagenames):
    """
        Creates pixel coordinate list from a physical coordinate
        list and a list of images.

        Parameters
        ---------
        coords
            A list of stacker.Coord coordinates.
        imagenames
            A list of images' paths.
    """

    pixcoords = CoordList(imagenames, 'pixel', unit='pix')

    for (i, imagename) in enumerate(imagenames):
        for coord in _getPixelCoords1Im(coords, imagename, i):
            pixcoords.append(coord)
    return pixcoords


def _getPixelCoords1Im(coords, imagename, imageNumber):
    from interval import interval
    import math

    try:
        from taskinit import ia
        ia.open(imagename)
        cs = ia.coordsys()
        Nx = ia.shape()[0]
        Ny = ia.shape()[1]
        ia.done()
        x0 = cs.referencevalue()['numeric'][0]
        y0 = cs.referencevalue()['numeric'][1]
        x_pix_ref = cs.referencepixel()['numeric'][0]
        y_pix_ref = cs.referencepixel()['numeric'][1]
        x_pix_inc = cs.increment()['numeric'][0]
        y_pix_inc = cs.increment()['numeric'][1]

# If we fail to load ia, we will use pyrap instead.
# This probably means stacker was loaded from outside casapy.
    except ImportError:
        from pyrap.images import image
        im = image(imagename)
        cs = im.coordinates().get_coordinate('direction')
        dir_axis_index = im.coordinates().get_axes().index(cs.get_axes())
        imshape = im.shape()
        try:
            x_axis_index = cs.get_axes().index('Right Ascension')
        except ValueError:
            raise ValueError('Could not find direction coordinate: '\
                              'RightAscension')
        try:
            y_axis_index = cs.get_axes().index('Declination')
        except ValueError:
            raise ValueError('Could not find direction coordinate: '\
                              'Declination')
        Nx = im.shape()[dir_axis_index+x_axis_index]
        Ny = im.shape()[dir_axis_index+y_axis_index]
        x0 = cs.get_referencevalue()[x_axis_index]
        y0 = cs.get_referencevalue()[y_axis_index]
        x_pix_ref = cs.get_referencepixel()[x_axis_index]
        y_pix_ref = cs.get_referencepixel()[y_axis_index]
        x_pix_inc = cs.get_increment()[x_axis_index]
        y_pix_inc = cs.get_increment()[y_axis_index]

    pixcoords = []

    for (cc,coord) in enumerate(coords):
        if imageNumber==coord.image:

            dx = (coord.x - x0)*math.cos(coord.y)
            dy = math.asin(math.sin(coord.y)/math.cos(dx)) - y0
            x = dx/x_pix_inc+x_pix_ref
            y = dy/y_pix_inc+y_pix_ref
            if x in interval[0., Nx-1] and y in interval[0., Ny-1]:
                c = Coord(x, y, z=coord.z, weight=coord.weight, image=imageNumber)
                try:
                    c.obsSpecArg = coord.obsSpecArg
                except AttributeError:
                    pass
                pixcoords.append(c)
            else:
                if x in interval[0, Nx-1]:
                    print x,y , Nx, Ny, coord, imagename, imageNumber
                    raise Exception(' y not in autorized interval')
                elif y in interval[0., Ny-1]:
                    print x,y , Nx, Ny, coord, imagename, imageNumber
                    raise Exception('x  not in autorized interval')
                else:
                    print x,y , Nx, Ny, coord, imagename, imageNumber
                    raise Exception('x and y not in autorized interval')
    return pixcoords
