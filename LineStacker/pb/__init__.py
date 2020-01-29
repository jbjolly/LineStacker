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

from __future__ import division
import math
import numpy as np
import stacker
import os

def guesspb(vis):
    """
    Produces a PrimaryBeamModel from a measurementset
    """
    print(vis)
    from taskinit import ms,tb,qa
    ms.open(vis)
    freq = (np.mean(ms.range('chan_freq')['chan_freq'])/1e9)
    ms.done()
    tb.open(vis+'/OBSERVATION')
    telescope = tb.getcol('TELESCOPE_NAME')[0]
    tb.done()
    pbfile = '{0}-{1:.1f}GHz.pb'.format(telescope, freq)
    print(pbfile)
    if not os.access(pbfile, os.F_OK):
        stacker.make_pbfile(vis, pbfile)
    return MSPrimaryBeamModel(pbfile)


class PrimaryBeamModel(object):
    """
    Base class for all primary beam models.
    Can be instanciated to get a constant beam model,
    i.e., set to 1 for all coordinates.
    """
    def __init__(self, maxdist=None, nu0=None, *args, **kwargs):
        """
        Constructor

        Keyword arguments:
        maxdist -- Maximum distance, pb(dx,dy) = 0 if x^2+y^2 > maxdist^2.
        nu0 -- Allows for frequency scaling of dx,dy -> (dx*nu0/nu, dy*nu0/nu)
        """

        super(PrimaryBeamModel, self).__init__(*args, **kwargs)
        self.maxdist = maxdist
        self.nu0 = nu0
        

    def __call__(self, dx, dy, nu=None):
        """
        Returns primary beam correction for translation (dx, dy) from center.

        Keyword arguments:
        dx -- Projected separation in x direction in radians
        dy -- Projected separation in y direction in radians
        nu -- Frequency, will apply no frequency scaling if set to None.
        """

        nuscaling = self._getnuscaling(nu)
        if self.maxdist is not None and (dx**2+dy**2)*nuscaling**2 > self.maxdist**2:
            return 0.

        return 1.


    def cdata(self):
        """
        Returns a version of the primary beam model that can be sent to c functions.
        """
        pbtype = stacker.PB_CONST
        pbfile = ''
        pbnpars = 0
        pbpars = None
        return pbtype, pbfile, pbnpars, pbpars


    def _getnuscaling(self, nu):
        nuscaling = 1.

        if self.nu0 is not None and nu is not None:
            nuscaling = self.nu0/nu

        return nuscaling


class GaussianPrimaryBeamModel(PrimaryBeamModel):
    """ 
    Primary beam model based on dish diameter
    """
    def __init__(self, dishdia='12m', nu0=None, *args, **kwargs):
        """
        Constructor

        Keyword arguments:
        dishdia -- Diamater of the telescope dish
        """
        super(GaussianPrimaryBeamModel, self).__init__(*args, **kwargs)
        self.dishdia = dishdia
        self.nu0 = nu0

    def __call__(self, dx, dy, nu=None):
        """
        Returns primary beam correction for translation (dx, dy) from center.

        Keyword arguments:
        dx -- Projected separation in x direction in radians
        dy -- Projected separation in y direction in radians
        nu -- Frequency, either self.nu0 or nu must be set.
        """
        from scipy.constants import c
        from taskinit import qa

        if nu is None:
            if self.nu0 is None:
                return 0.
            nu = self.nu0

        dishdia = qa.convert(self.dishdia, 'm')['value']
        vp_fwhm = 1.22*c/nu/dishdia
        return np.exp(-4.*np.log(2)*(dx**2+dy**2)/vp_fwhm**2)


class MSPrimaryBeamModel(PrimaryBeamModel):
    """ 
    Primary beam model based on casa image. 
    """

    def __init__(self, imagename, *args, **kwargs):
        """
        Constructor

        Keyword arguments:
        imagename -- Str to casa image of primary beam
        """
        super(MSPrimaryBeamModel, self).__init__(*args, **kwargs)

        self.imagename = imagename

        try:
            from taskinit import ia
            ia.open(imagename)
            self.cs = ia.coordsys()

            self.nx = ia.shape()[0]
            self.ny = ia.shape()[1]
            self.refpix_x    = self.cs.referencepixel()['numeric'][0]
            self.refpix_y    = self.cs.referencepixel()['numeric'][1]
            self.increment_x = self.cs.increment()['numeric'][0]
            self.increment_y = self.cs.increment()['numeric'][1]

            try:
                self.frequencyaxis = self.cs.findaxisbyname('frequency')
                self.nu0 = self.cs.referencevalue()['numeric'][self.frequencyaxis]
            except Exception:
                self.nu0 = None
                print('Some stuff!')

            self.data = ia.getregion()[:,:,0,0]
            ia.done()
        except ImportError:
            from pyrap.images import image
            im = image(imagename)
            self.nx = im.shape()[-1]
            self.nx = im.shape()[-2]
            self.cs = im.coordinates()
            self.cs_dir = self.cs.get_coordinate('direction')

            self.refpix_x = self.cs_dir.get_referencepixel()[1]
            self.refpix_y = self.cs_dir.get_referencepixel()[0]
            self.increment_x = self.cs_dir.get_increment()[1]
            self.increment_y = self.cs_dir.get_increment()[0]
            try:
                self.nu0 = self.cs.get_coordinate('spectral').get_referencevalue()
            except Exception:
                self.nu0 = None
                print('Warning! No frequency information in primary beam model.')

            self.data = im.getdata()[0,0]


    def __call__(self, dx, dy, nu = None):
        """
        Returns primary beam correction for translation (dx, dy) from center.

        Keyword arguments:
        dx -- Projected separation in x direction in radians
        dy -- Projected separation in y direction in radians
        nu -- Frequency, will apply no frequency scaling if set to None.
        """
        nuscaling = self._getnuscaling(nu)
        x = int(nuscaling*(dx)/self.increment_x + self.refpix_x + 0.5)
        y = int(nuscaling*(dy)/self.increment_y + self.refpix_y + 0.5)

        if x < 0 or x >= self.nx:
            return 0.

        if y < 0 or y >= self.ny:
            return 0.

        if self.data is None:
            return 0.

        return self.data[x,y]


    def cdata(self):
        """
        Returns a version of the primary beam model that can be sent to c functions.
        """
        pbtype = stacker.PB_MS
        pbfile = self.imagename
        pbnpars = 0
        pbpars = None
        return pbtype, pbfile, pbnpars, pbpars

