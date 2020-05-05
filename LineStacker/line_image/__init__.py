"""
**LineStacker.line_image Module:**\n
Module for cube line stacking.
"""
# LineStacker, Python module for stacking of interferometric data.
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
import LineStacker
import LineStacker.image
import LineStacker.tools
import LineStacker.analysisTools
import numpy as np
from sys import modules


skymap = []
data = []
oldimagenames = []
stampsize = 0
imagesizes = []


def stack(  coords,
            outfile='stackResult.image',
            stampsize = 32,
            imagenames= [],
            method = 'mean',
            spectralMethod='z',
            weighting = None,
            maskradius=0,
            psfmode = 'point',
            primarybeam = None,
            fEm = 0,
            chanwidth='default',
            plotIt=False,
            regridFromZ=False,
            regridMethod='scaleToMin',
            saveSubCubes=False,
            **kwargs):
    """
        Performs line stacking in the image domain.
        returns: Estimate of stacked flux assuming point source.

        Parameters
        ---------
        coords
            A coordList object of all target coordinates.
	    outfile
            Target name for stacked image.
        stampsize
            size of target image in pixels
        imagenames
            Name of images to extract flux from.
        method
            'mean' or 'median', will determined how pixels are calculated
        spectralMethod
            Method to select the central frequency of the stack. The corresponding value should be found in the 3rd column of the coord file. Possible methods are:\n
            **'z'**: the redshift of the observed line, if used
            **fEm** (emission frequency) must be informed as an argument of this Stack function.\n
            **'centralFreq'**: the (observed) central frequency, or velocity.\n
            **'channel'**: to dirrectly input the channel number of the center of the stack.
        weighting
            used if method set to **'mean'**, possible values are:\n
            **'sigma2'**, weights are set to 1/sigma**2 where sigma is the standard deviation of the corresponding data cube (excluding masked region). See **calculate_sigma2_weights**\n
            **'sigma2F'**, similar to **sigma2** except weights are individually computed for each spectral bin. See **calculate_sigma2_weights_spectral**\n
            **'1/A'**, **/!\\\\ only use if lines are visible pre-stacking**, weights are set to one over the line amplitude. See **calculate_amp_weights**\n
            **None**, using weights in coords (set to 1 by default).
        maskradius
            allows blanking of centre pixels in weight calculation
        psfmode
            Allows application of filters to stacking, currently not supported.
        primarybeam
            only applies if weighting='pb'
        fEm
            rest emission frequency of the line,
        chanwidth
            number of channels of the resulting stack, default is number of channels in the first image.
        plotIt
            direct plot option
        regridFromZ
            if set to True new images will be created, regridded to take into account the redshift difference of the different sources. See LineStacker.analysisTools.regridFromZ for a more complete description.
            NB: ALL IMAGES SHOULD HAVE SAME FREQUENCY BIN originally.
        regridMethod
            Used if regridFromZ is True. Can be set either to 'scaleToMin' or 'scaleToMax'.
            In the first case all images are regrided to match the smallest redshift (over-gridding), all images are regridded to match the highest redshift in the other case (under-gridding)
        saveSubCubes
            If set to True (or to str) sub cubes will be saved as a numpy file. (outSubCubes.npy if set to True, user defined is set to str)
    """

    from ..interval import interval
    import os
    import shutil
    import numpy as np
    from taskinit import ia, casalog

    casalog.origin('line-stacker')
    casalog.post('#'*42,'INFO')
    casalog.post('#'*5 +  ' {0: <31}'.format("Begin Task: Line-Stacker")+'#'*5, 'INFO')

    global skymap
    global data
    global oldimagenames


    if coords.coord_type == 'physical':
        coords = LineStacker.getPixelCoords(coords, imagenames)

    if chanwidth=='default':
        ia.open(imagenames[0])
        chanwidth=ia.summary()['shape'][-1]
        ia.done()

    if regridFromZ:

        ia.open(imagenames[0])
        imInfo=ia.summary()
        ia.done()

        if imInfo['axisnames'][-1]=='Frequency':
            imagenames, coords=LineStacker.analysisTools.regridFromZ(
                                                    coords,
                                                    imagenames,
                                                    stampsize=stampsize,
                                                    chanwidth=chanwidth,
                                                    fEm=fEm,
                                                    regridMethod=regridMethod)
        else:
            raise Exception('last axis is not Frequency, cant regrid')
    #fills skymap and loads data as zeros
    _allocate_buffers(coords.imagenames, stampsize, len(coords), chanwidth)

    #fills data with skymap accordingly
    _load_stack(coords, psfmode, fEm=fEm, spectralMethod=spectralMethod)
    if method=='mean' and weighting=='sigma2':
        calculate_sigma2_weights(coords, maskradius=maskradius)
    if method=='mean' and weighting=='sigma2F':
        calculate_sigma2_weights_spectral(coords, maskradius=maskradius)
    if method=='mean' and weighting=='1/A':
        calculate_amp_weights(coords, **kwargs)

    #actual stack
    stacked_im  = _stack_stack(method, coords)


    #write image output
    _write_stacked_image(outfile, stacked_im,
                         coords,
                         stampsize,
                         chanwidth,
                         fEm,
                         regridFromZ=regridFromZ,
                         regridMethod=regridMethod)
    if saveSubCubes!=False:
        if type(saveSubCubes)==str:
            np.save(saveSubCubes+'.npy', data)
            outCoordFile=open(saveSubCubes+'_outCoordsFile.txt','w')
            for (co, coord) in enumerate(coords):
                outCoordFile.write(str(coord.x)+', '+str(coord.y)+', '+str(coord.z)+'\n')
            outCoordFile.close()
        else:
            np.save('outSubCubes.npy', data)
            outCoordFile=open('outCoordsFile.txt','w')
            for (co, coord) in enumerate(coords):
                outCoordFile.write(str(coord.x)+', '+str(coord.y)+', '+str(coord.z)+'\n')
            outCoordFile.close()
    #plt.figure()
    if plotIt:
        import matplotlib.pyplot as plt
        try:
            stackedImIncludingBorderEffects=stacked_im[int(stampsize/2), int(stampsize/2),0,:]*(len(coords)/(len(coords)-unusedFrequencies))
        except ZeroDivisionError:
            print 'ERORR : your chanwidth is too large and no frequencies were found for ANY of your images in one or more bins'

        plt.figure()
        #ax=im.add_subplot(2)
        plt.plot(stacked_im[int(stampsize/2), int(stampsize/2),0,:], 'r', linewidth=5)
        plt.plot(unusedFrequencies*max(stacked_im[int(stampsize/2), int(stampsize/2),0,:])/len(coords), 'b', linewidth=2.5, label='percentage of the sources not used')
        plt.plot(stackedImIncludingBorderEffects, 'g', linewidth=5, label='corrected')
        plt.tick_params(axis='x', labelsize=25)
        plt.tick_params(axis='y', labelsize=20)
        plt.xlabel('Rest Frequency (GHz)', size=40)
        #plt.ylabel('Flux density (Jy)', size=25)
        plt.ticklabel_format(style='plain')
        plt.title('Stacked line', size=40)
        plt.yticks(np.arange(-0.0001, 0.0006, 0.0002))
        plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
        plt.legend()
        plt.show()
    #
    casalog.post('#'*5 +  ' {0: <31}'.format("End Task: stacker")+'#'*5)
    casalog.post('#'*42)
    oldimagenames=imagenames
    return [stacked_im[int(stampsize/2), int(stampsize/2),0,int(chanwidth/2)], stacked_im[int(stampsize/2), int(stampsize/2),0,:], stacked_im, unusedFrequencies]


def _allocate_buffers(  imagenames,
                        new_stampsize,
                        nstackpos,
                        new_chanwidth):

    #buffers skymap (all cubes) and data (empty cubes that will be filled and stacked)
    try:
        from taskinit import ia
    except ImportError:
        raise Exception('ia import failed, please run Line-Stacker from CASA')

    global skymap
    global data
    global oldimagenames
    global stampsize
    global imagesInfo
    global chanwidth

    #ia.open(imagenames[0])
    #cs = ia.coordsys()
    #outnchans = ia.boundingbox()['trc'][3]+1
    #ia.done()

    outnstokes = 1
    chanwidth=new_chanwidth
# To improve performance this module will keep buffers between run.
# This following code resets these buffers if they have grown obsolete.
    if oldimagenames == []:
            oldimagenames = imagenames

    if oldimagenames != imagenames:
            oldimagenames = imagenames
            skymap = []
            data = []

    if stampsize == 0:
            stampsize = new_stampsize

    elif stampsize != new_stampsize:
            stampsize = new_stampsize
            data = []
            skymap = []

    if not(data == []) and nstackpos != data.shape[0]:
        data = []

# If there is no skymap buffer create one.
# This is the data that is most important to buffer.
# Reading a skymap from disk is a time consuming task and we don't want to do this too much.
    if skymap == []:
        imagesInfo={}
        imagesInfo['imagesizes']=[]
        imagesInfo['freq0']=[]
        imagesInfo['freqBin']=[]
        imagesInfo['freqlen']=[]
        for (i,imagename) in enumerate(imagenames):

            ia.open(imagename)

            skymap.append(ia.getregion())
            infos = ia.summary()
            imagesInfo['imagesizes'].append((ia.shape()[0], ia.shape()[1],ia.shape()[3]))
            imagesInfo['freq0'].append(infos['refval'][3]-infos['refpix'][3]*infos['incr'][3])
            imagesInfo['freqBin'].append(infos['incr'][3])
            imagesInfo['freqlen'].append(ia.shape()[3])

            ia.done()
    else:
        for imagename in imagenames:
            ia.open(imagename)
            infos = ia.summary()
            imagesInfo['freq0'].append(infos['refval'][3])
            imagesInfo['freqBin'].append(infos['incr'][3])
            imagesInfo['imagesizes'].append((ia.shape()[0], ia.shape()[1],ia.shape()[3]))
            imagesInfo['freqlen'].append(ia.shape()[3])
            ia.done()

    # If there is no data buffer create one.
    # The data buffer is used to save the right stacking positions before stacking them.
    # During stacking this is where the full stack will actually be saved.
    if data == []:
        data = np.zeros((nstackpos, new_stampsize, new_stampsize, outnstokes, new_chanwidth))
    else:
        data = 0.*data

def calculate_sigma2_weights(coords, maskradius=0.):
    """
        Computes standard deviation of data cubes and sets weights to one over sigma**2.

        Parameters
        ---------
        coords
            A coordList object of all target coordinates.
        maskradius
            Radius (in pixel) of the mask, centered on coordinate center, to avoid including pixels close to source in noise computation.
    """
    if maskradius>=stampsize:
        raise Exception('the mask radius is bigger than the stamp size')
    if stampsize<=4:
        print 'stamp size is very small, using sigma2 weight is probably irrelevant'

    if maskradius!=0:
        X = np.arange(0, stampsize)-int(stampsize/2)
        Y = np.arange(0, stampsize)-int(stampsize/2)
        X,Y = np.meshgrid(X,Y)

    for i,coord in enumerate(coords):
        tmpdata = np.copy(data[i,:,:,:,:])
        if maskradius!=0:
            for j in range(tmpdata.shape[2]):
                for k in range(tmpdata.shape[3]):
                    tmpdata[:,:,j,k]  = (tmpdata[:,:,j,k]*np.double( X**2+Y**2>maskradius**2))
        else:
            pass
        sigma = np.std(tmpdata[np.nonzero(tmpdata)])
        if sigma == 0:
            coord.weight = coord.weight*0.
        else:
            coord.weight = coord.weight*1./sigma**2

def calculate_sigma2_weights_spectral(coords, maskradius=0.):
    """
        Computes standard deviation of data cubes in every spectral channel and sets weights to one over sigma**2.

        Parameters
        ---------
        coords
            A coordList object of all target coordinates.
        maskradius
            Radius (in pixel) of the mask, centered on coordinate center, to avoid including pixels close to source in noise computation.
    """

    if maskradius>=stampsize:
        raise Exception('the mask radius is bigger than the stamp size')
    if stampsize<=4:
        print 'stamp size is very small, using sigma2 weights may be biased'
    if maskradius!=0:
        X = np.arange(0, stampsize)-int(stampsize/2)
        Y = np.arange(0, stampsize)-int(stampsize/2)
        X,Y = np.meshgrid(X,Y)
    for i,coord in enumerate(coords):
        tmpdata = np.copy(data[i,:,:,:,:])
        sigmaF=np.zeros(tmpdata.shape[3])
        if maskradius!=0:
            for k in range(tmpdata.shape[3]):
                for j in range(tmpdata.shape[2]):
                    tmpdata[:,:,j,k]  = (tmpdata[:,:,j,k]*np.double( X**2+Y**2>maskradius**2))
                sigmaF[k]=np.std(tmpdata[:,:,0,k])#[np.nonzero(tmpdata[:,:,0,k])])
        else:
            for k in range(tmpdata.shape[3]):
                sigmaF[k]=np.std(tmpdata[:,:,0,k])
        try:
            for k in range(len(coord.weight)):
                if coord.weight[k]==0:
                    pass
                else:
                    coord.weight[k]=1./(sigmaF[k]*sigmaF[k])

        except TypeError:
            coord.weight=1./(sigmaF*sigmaF)

def calculate_amp_weights(coords, fit=False):
    """
        Sets weights of each coord to one over the line amplitude.\n
        **/!\\\\ only use if lines are visible pre-stacking**

        Parameters
        ---------
        coords
            A coordList object of all target coordinates.
        fit
            If set to **True**, spectra will be extracted from each coordinate pixel and fited with a gaussian. The gaussian's extracted amplitude will be used as the line's amplitude.\n
            If set to **False**, the line's amplitude is set to the value of the brightest spectral bin of each coordinate pixel.
    """

    if not fit:
        import numpy as np
        for i,coord in enumerate(coords):
            amp=np.max(data[i,int(data.shape[1]/2.), int(data.shape[2]/2.), 0, :])
            coord.weight=coord.weight*1./amp
    else:
        import LineStacker.tools.fit as fitTool
        for i,coord in enumerate(coords):
            fitos=fitTool.GaussFit(fctToFit=data[i,int(data.shape[1]/2.), int(data.shape[2]/2.), 0, :], returnInfos=True)
            amp=fitos[1][0]
            coord.weight=coord.weight*1./amp

def _load_stack(coords, psfmode='point', fEm=0, spectralMethod='z', Return=False):
    from ..interval import interval

    global data
    global unusedFrequencies
    listWeights=False
    unusedFrequencies=np.zeros(chanwidth)
    if len(coords) > data.shape[0]:
        _allocate_buffers(coords.imagenames, stampsize, len(coords))
    #the number of cubes not used at the given spectral bin
    for (i,coord) in enumerate(coords):
        if spectralMethod=='z':
            obsFreq=fEm/(1.+coord.z)
            obsFreqArg=int(round((obsFreq-imagesInfo['freq0'][coord.image])/imagesInfo['freqBin'][coord.image]))
        elif spectralMethod=='centralFreq':
            obsFreq=coord.z
            obsFreqArg=int(round((obsFreq-imagesInfo['freq0'][coord.image])/imagesInfo['freqBin'][coord.image]))
        elif spectralMethod=='channel':
            obsFreqArg=int(round(coord.z))
        #coord.obsSpecArg=obsFreqArg
        if coord.obsSpecArg!=0:
            obsFreqArg=coord.obsSpecArg
        blcx=int(coord.x - int(stampsize/2.)+0.5)
        blcy=int(coord.y - int(stampsize/2.)+0.5)
        blcf = obsFreqArg - int(chanwidth/2.)
        numberOfZeroLeftF=0

        trcx = blcx + stampsize
        trcy = blcy + stampsize
        trcf = blcf + chanwidth
        numberOfZeroRightF=0
        #coords[i].setWeightToX(chanwidth,X=1)
        if blcf<0:
            listWeights=True
            numberOfZeroLeftF=blcf
            coords[i].setZeroWeightLeft(numberOfZeroLeftF,imagesInfo['freqlen'][coord.image])
            if len(coord.weight)>chanwidth:
                coord.weight=coord.weight[:chanwidth]
            blcf=0

        if trcf>imagesInfo['freqlen'][coord.image]-1:
            listWeights=True
            numberOfZeroRightF=trcf-imagesInfo['freqlen'][coord.image]
            coords[i].setZeroWeightRight(numberOfZeroRightF, imagesInfo['freqlen'][coord.image])
            trcf=imagesInfo['freqlen'][coord.image]
            if len(coord.weight)>chanwidth:
                coord.weight=coord.weight[len(coord.weight)-chanwidth:]


# Currently including positions on the edge of the skymap (in space) will raise an error.
# This could certainly be improved upon.
        if (interval[blcx, trcx] in interval[0, imagesInfo['imagesizes'][coord.image][0]]
                and interval[blcy, trcy] in interval[0, imagesInfo['imagesizes'][coord.image][1]]):

            if not numberOfZeroRightF and not numberOfZeroLeftF:
                data[i,:,:,0,:] = skymap[coord.image][blcx:trcx,blcy:trcy,0,blcf:trcf]
            elif numberOfZeroRightF:
                data[i,:,:,0,-numberOfZeroLeftF:-numberOfZeroRightF] = skymap[coord.image][blcx:trcx,blcy:trcy,0,blcf:trcf]
                unusedFrequencies[0:-numberOfZeroLeftF]+=1
                unusedFrequencies[-numberOfZeroRightF:]+=1

            elif numberOfZeroLeftF and not numberOfZeroRightF:
                data[i,:,:,0,-numberOfZeroLeftF:] = skymap[coord.image][blcx:trcx,blcy:trcy,0,blcf:trcf]
                unusedFrequencies[0:-numberOfZeroLeftF]+=1

        else:
            raise Exception('the source located at '+str(coord.x)+' '+str(coord.y)+' on image '+str(coord.image)+' is too close to the edge, trying to stack out of boundaries, try stacking with smaller stampsize')

            data[i,:,:,0,:] = 0
        if psfmode == 'star':
            raise Exception('star not supported')
    #fix so that all weights are list like if some weights have been turned into list due to zeros left of right
    if listWeights:
        for coord in coords:
            try:
                len(coord.weight)
            except TypeError:
                coord.weight=np.array([coord.weight for i in range(chanwidth)])
    if Return:
        return data

def _write_stacked_image(   imagename,
                            pixels,
                            coords,
                            stampsize,
                            chanwidth,
                            fEm,
                            regridFromZ=False,
                            regridMethod='scaleToMin'):
    import os
    import shutil
    from taskinit import ia
#     global stampsize
    if os.access(imagename, os.F_OK): shutil.rmtree(imagename)


    template_image=coords.imagenames[0]
    ia.open(template_image)
    beam = ia.restoringbeam()
    cs = ia.coordsys()
    ia.done()

    csnew = cs.copy()
    csnew.setreferencevalue([0.]*2, 'dir')
    csnew.setreferencepixel([int(stampsize/2+0.5)]*2, 'dir')
    if not regridFromZ:
        if fEm==0:
            centralFreq=cs.referencevalue()['numeric'][3]
        else:
            centralFreq=fEm
        csnew.setreferencevalue(centralFreq, type='spectral')

    elif regridFromZ:
        c=299792.458#in km
    	if regridMethod=='scaleToMin':
            nonZeroZ=[]
            nonZeroIndex=[]
            for (i, coord) in enumerate(coords):
                if coord.z==0:
                    pass
                else:
                    nonZeroZ.append(coord.z)
                    nonZeroIndex.append(i)
            zToScale=min([zz for zz in nonZeroZ])
            scaleIndex=np.argmin([zz for zz in nonZeroZ])
            scaleIndex=nonZeroIndex[scaleIndex]
        elif regridMethod=='scaleToMax':
            zToScale=max([cc.z for cc in coords])
            scaleIndex=np.argmax([cc.z for cc in coords])

        csnew.setreferencevalue(float(fEm)/(1+zToScale), type='spectral')
        ia.open(coords.imagenames[scaleIndex])
        newVelIncr=ia.summary()['incr'][-1]
        ia.done()
        csnew.setincrement(value=newVelIncr, type='spectral')

    csnew.setreferencepixel([int(chanwidth/2+0.5)]*2, type='spectral')
    ia.fromarray(imagename, pixels=pixels, csys = csnew.torecord())
    ia.open(imagename)
    #ia.setrestoringbeam(beam=beam)
    ia.done()

def  Weights_freq_image(imagename):
    from taskinit import ia, casalog

    ia.open(imagename)
    cs1=ia.coordsys()
    allPix=ia.getchunk()
    ia.done()
    thisImWeights=np.zeros(allPix.shape[-1])
    for freqbin in range(allPix.shape[-1]):
        sigma=np.std(allPix[:,:,:,freqbin])
        thisImWeights[freqbin]=1./(sigma*sigma)
    return thisImWeights

def myVeryOwnWeightedAverage(data,coords):

    #if weights for each frequencies and not just for each image
    allPix=np.zeros(data.shape[1:])
    for i in range(chanwidth):
        freqWeight=([coord.weight[i] for coord in coords])
        #if freqWeight==([0 for coord in coords]) or freqWeight==([0.0 for coord in coords]) or freqWeight is np.zeros(len(coords)):
        #if np.array(freqWeight) is np.zeros(len(coords)):
        #    tempPix=0
        #else:
        if np.sum(freqWeight)!=0:
            tempPix=np.average(data[0:len(coords),:,:,:,i], 0, weights=freqWeight)
        else:
            tempPix=np.zeros(allPix.shape[:-1])
        #allPix.append(tempPix)
        allPix[:,:,:,i]=tempPix
    return allPix


def _stack_stack(method, coords):
    """
        Performs the actual stacking on the data in the stack.
        All data should be loaded in to stack before calling this function.
    """
    #pixles array will be filled with the stack values
    import matplotlib.pyplot as plt
    pixels = np.zeros(data.shape[1:])
    if method == 'median':
        pixels = np.median(data[0:len(coords),:,:,:,:], 0)

    elif method == 'mean':
        try:
            len(coords[0].weight)
            pixels=myVeryOwnWeightedAverage(data,coords)
        except TypeError:
            pixels = np.average(data[0:len(coords),:,:,:,:], axis=0, weights=([coord.weight for coord in coords]))
    return pixels
