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
# Boston, MA  02110-1301, USA.
from ctypes import c_double, POINTER, c_char_p, c_int, c_bool
import numpy as np
import stacker
import stacker.pb

c_stack = stacker.libstacker.stack
c_stack.restype = c_double
c_stack.argtype = [c_int, c_char_p, c_int,
                   c_int, c_char_p, c_int,
                   c_int, c_char_p, POINTER(c_double), c_int,
                   POINTER(c_double), POINTER(c_double), POINTER(c_double),
                   c_int, c_bool]
c_stack_mc = stacker.libstacker.stack_mc
c_stack_mc.argtype = [c_int, c_char_p, c_int,
                      c_int, c_char_p, POINTER(c_double), c_int,
                      POINTER(c_double), POINTER(c_double), POINTER(c_double),
                      c_int, c_int, POINTER(c_char_p),
                      POINTER(c_double), POINTER(c_double), c_int, c_bool]


def stack(coords, vis, outvis='', imagename='', cell='1arcsec', stampsize=32,
          primarybeam='guess', datacolumn='corrected', use_cuda = False):
    """
         Performs stacking in the uv domain.


         coords      -- A coordList object of all target coordinates.
         vis         -- Input uv data file.
         outvis      -- Output uv data file. Can be set to '' to not save
                        stacked visibilities.
         datacolumn  -- Either 'corrected' or 'data'. Which column stacking is
                        applied to.
         primarybeam -- How to calculated primary beam. Currently only two
                        options, 'guess' (using casa builtin model) or
                        'constant' (i.e. no correction)
         imagename   -- Optional argument to image stacked data.
         cell        -- pixel size for target image
         stampsize   -- size of target image in pixels

         returns: Estimate of stacked flux assuming point source.
    """
    import shutil
    import os
    try:
        from taskinit import casalog
    except ImportError:
        casalog = None

    if casalog is not None:
        casalog.origin('stacker')
        casalog.post('#'*42, 'INFO')
        casalog.post('#'*5 + ' {0: <31}'.format("Begin Task: Stacker")+'#'*5,
                     'INFO')
        casalog.post('Number of stacking positions: {0}'.format(len(coords)),
                     'INFO')

    if outvis != '':
        if not os.access(outvis, os.F_OK):
            shutil.copytree(vis, outvis)

    infiletype, infilename, infileoptions = stacker._checkfile(vis, datacolumn)
    if outvis != '':
        outfiletype, outfilename, outfileoptions =\
            stacker._checkfile(outvis, datacolumn)
    else:
        outfilename = ''
        outfiletype = stacker.FILE_TYPE_NONE
        outfileoptions = 0

    if casalog is not None:
        casalog.post('Input uv file: \'{0}\' of type {1}'.format(
            infilename, stacker.FILETYPENAME[infiletype]), 'INFO')
        if outvis != '':
            casalog.post('Output uv file: \'{0}\' of type {1}'.format(
                outfilename, stacker.FILETYPENAME[outfiletype]), 'INFO')
        else:
            _ = 'No output uv file given, will not write stacked visibility'
            casalog.post(_, 'INFO')

# primary beam
    if primarybeam == 'guess':
        primarybeam = stacker.pb.guesspb(vis)
    elif primarybeam in ['constant', 'none'] or primarybeam is None:
        primarybeam = stacker.pb.PrimaryBeamModel()

    pbtype, pbfile, pbnpars, pbpars = primarybeam.cdata()

    x = [p.x for p in coords]
    y = [p.y for p in coords]
    weight = [p.weight for p in coords]

    x = (c_double*len(x))(*x)
    y = (c_double*len(y))(*y)
    weight = (c_double*len(weight))(*weight)

    import time
    start = time.time()
    flux = c_stack(infiletype, c_char_p(infilename), infileoptions,
                   outfiletype, c_char_p(outfilename), outfileoptions,
                   pbtype, c_char_p(pbfile), pbpars, pbnpars,
                   x, y, weight, c_int(len(coords)), c_bool(use_cuda))
    stop = time.time()
#     print("Started stack at {}".format(start))
#     print("Finished stack at {}".format(stop))
    print("Time used to stack: {0}".format(stop-start))

    if imagename != '':
        import clean
        import clearcal
        clearcal.clearcal(vis=outvis)
        clean.clean(vis=outvis, imagename=imagename, field='0', mode='mfs',
                    cell=cell, imsize=stampsize, weighting='natural')

    if casalog is not None:
        casalog.post('#'*5 + ' {0: <31}'.format("End Task: stacker")+'#'*5)
        casalog.post('#'*42)
    return flux


def noise(coords, vis, weighting='sigma2', imagenames=[], beam=None,
            nrand=50, stampsize=32, maskradius=None, datacolumn='corrected'):
    """ Calculate noise using a Monte Carlo method, can be time consuming. """
    import stacker
    import stacker.image
    from math import pi
    print 'a1'
    if beam is None:
        try:
            from taskinit import ia, qa

            ia.open(imagenames[0])
            beam = qa.convert(ia.restoringbeam()['major'], 'rad')['value']
            ia.done()
        except ImportError:
            beam = 1/3600./180.*pi
    print 'a2'
    dist = []
    '''
    /1\ changer dist=[] en dist=np.zeros(nrand)
    '''
    for i in range(nrand):
        random_coords = stacker.randomizeCoords(coords, beam=beam)
        print i, 'coords rand'
        if weighting == 'sigma2':
            print 'sig2'
            random_coords = stacker.image.calculate_sigma2_weights(
                random_coords, imagenames, stampsize, maskradius)
        laStackence=stack(random_coords, vis, datacolumn=datacolumn)
        #dist.append(stack(random_coords, vis))
        dist.append(laStackence)

        print 'append donc'
    '''
    /1\ =v ne faudrait il pas return egalement le mean ? =v
    '''
    print 'a3'
    np.savetxt('noise.txt',np.real(np.array(dist)))
    return [np.std(abs(np.real(np.array(dist)))), np.mean(abs(np.real(np.array(dist))))]


def noise_fast(coords, models, vis, datacolumn='corrected',
               primarybeam='guess', use_cuda=True,
               nbin=None, bins=None):
    """ Calculate noise using a Monte Carlo method, can be time consuming. """
    import stacker

    if not use_cuda:
        raise NotImplementedError
    if len(coords) != len(models):
        raise RuntimeError('Number of coordinate objects does not match number of models.')
    nmc = len(coords)
    infiletype, infilename, infileoptions = stacker._checkfile(vis, datacolumn)

    if primarybeam == 'guess':
        primarybeam = stacker.pb.guesspb(vis)
    elif primarybeam in ['constant', 'none'] or primarybeam is None:
        primarybeam = stacker.pb.PrimaryBeamModel()
    pbtype, pbfile, pbnpars, pbpars = primarybeam.cdata()

    coordlistlen = len(coords[0])
    for coordlist in coords:
        if coordlistlen != len(coordlist):
            raise RuntimeError('Number of coordinates must be same in all samples.')
    x = []
    y = []
    weight = []
    for coordlist in coords:
        x.extend([p.x for p in coordlist])
        y.extend([p.y for p in coordlist])
        weight.extend([p.weight for p in coordlist])
    print('len(x) = {}'.format(len(x)))

    x = (c_double*len(x))(*x)
    y = (c_double*len(y))(*y)
    weight = (c_double*len(weight))(*weight)

    c_models = (c_char_p*len(models))(*models)

    if bins is None or nbin is None:
        raise 'No bins given!'
    if len(bins) != nbin+1:
        raise 'Number of bins must match nbin!'
    c_bins = (c_double*(nbin+1))(*bins)
    res_flux   = (c_double*(nmc*nbin))(*([0]*(nmc*nbin)))
    res_weight = (c_double*(nmc*nbin))(*([0]*(nmc*nbin)))

    c_stack_mc(infiletype, c_char_p(infilename), infileoptions,
               pbtype, c_char_p(pbfile), pbpars, pbnpars,
               x, y, weight, c_int(len(coords[0])), c_int(nmc),
               c_models,
               res_flux, res_weight, c_bins, c_int(nbin),
               use_cuda)

    return np.array(list(res_flux)), np.array(list(res_weight))
