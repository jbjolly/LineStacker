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
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

from ctypes import c_double, c_bool, POINTER, c_char_p, cdll, c_int
import numpy as np
import stacker
import stacker.pb

# lib = cdll.LoadLibrary('libmodsub.so')
lib = stacker.libstacker
c_modsub = lib.modsub
c_modsub.restype = c_int
c_modsub.argtype = [c_char_p, c_char_p, c_char_p,c_char_p]
c_modsub.argtype = [c_int, c_char_p, c_int,
                    c_int, c_char_p, c_int,
                    c_char_p,
                    c_int, c_char_p, POINTER(c_double), c_int,
                    c_bool, c_bool]

def modsub(model, vis, outvis='', datacolumn='corrected', primarybeam='guess', subtract=True, use_cuda=False, field = None):
    import shutil
    import os

	
    if outvis != '':
        if not os.access(outvis, os.F_OK):
            shutil.copytree(vis, outvis)

    print 1		
# primary beam
    if primarybeam == 'guess':
        primarybeam = stacker.pb.guesspb(vis)
    elif primarybeam in ['constant', 'none'] or primarybeam is None:
        primarybeam = stacker.pb.PrimaryBeamModel()
    pbtype, pbfile, pbnpars, pbpars = primarybeam.cdata()
    print 2
    if field is not None:
        select_field = True
        field = str(field)
    else:
        select_field = False
        field = ''
    print 3
    infiletype, infilename, infileoptions = stacker._checkfile(vis, datacolumn)
    print 4	    
    outfiletype, outfilename, outfileoptions = stacker._checkfile(outvis, datacolumn)
    print 5
    flux = c_modsub(infiletype, c_char_p(infilename), infileoptions,
                    outfiletype, c_char_p(outfilename), outfileoptions,
                    c_char_p(model), 
                    pbtype, c_char_p(pbfile), pbpars, pbnpars,
                    c_bool(subtract), c_bool(use_cuda),
                    c_bool(select_field), c_char_p(field))
    return 0


def cl_from_im(image, clname=None, threshold=None):
    from taskinit import qa,ia, cl

    ia.open(image)
    data = ia.getregion()
    cs = ia.coordsys()
    ia.done()

    data = np.where(np.isnan(data), 0, data)
    if threshold is None:
        datanz = np.nonzero(data)
    else:
        datanz = np.nonzero(np.abs(data) > threshold)

    modelfluxes = data[datanz]
    modelpixels = np.array(datanz)
    modellist = cs.convertmany(modelpixels, 
                               unitsin=['pix', 'pix', 'pix', 'pix'], 
                               unitsout=['rad', 'rad', '', 'Hz'])

    cl.done()

    for i in range(modellist.shape[1]):
        x = qa.formxxx(str(modellist[0, i])+'rad', format='hms', prec=6)
        y = qa.formxxx(str(modellist[1, i])+'rad', format='dms', prec=6)
        pos = ' '.join(['J2000', x, y])
        freq = str(modellist[3, i])+'Hz'
        flux = modelfluxes[i]
        cl.addcomponent(flux=flux, fluxunit='Jy', dir=pos, freq=freq)

    if clname is not None:
        cl.rename(clname)
        cl.done()
    else:
        return cl

        
