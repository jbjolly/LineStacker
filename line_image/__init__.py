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
import numpy as np
from sys import modules


skymap = []
data = []
oldimagenames = []
stampsize = 0
imagesizes = []


def stack(  coords,
            outfile,
            stampsize = 32,
            imagenames= [],
            method = 'mean',
            spectralMethod='z',
            weighting = None,
            maskradius=0,
            psfmode = 'point',
            primarybeam = None,
            fEm = 0,
            chanwidth=30,
            plotIt=False):
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
            How to select the central frequency of the stack\n
            The corresponding value should be found in the 3rd column of the coord file\n
            possible methods are:\n
            'z': the redshift of the observed line,
            additionaly fEm (emission frequency) must be informed (as an argument of this Stack function)
            'centralFreq': the (observed) central frequency, or velocity\n
            'channel': to dirrectly input the channel number of the center of the stack
        weighting
            only for method 'mean', if set to None will use weights in coords.
        maskradius
            allows blanking of centre pixels in weight calculation
        psfmode
            Allows application of filters to stacking, currently not supported.
        primarybeam
            only applies if weighting='pb'
        fEm
            rest emission frequency of the line,
        chanwidth
            number of channels of the resulting stack,
        plotIt
            direct plot option):
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

    #fills skymap and loads data as zeros
    _allocate_buffers(coords.imagenames, stampsize, len(coords), chanwidth)

    if method=='mean' and coords[0].weight==1 and weighting=='sigma2':
        _calculate_sigma2_weights(coords, maskradius=maskradius)
    if method=='mean' and coords[0].weight==1 and weighting=='sigma2F':
        _calculate_sigma2_weights_spectral(coords, maskradius=maskradius)
    if method=='mean' and coords[0].weight==1 and weighting=='1/A':
        _calculate_amp_weights(coords)

    #fills data with skymap accordingly
    _load_stack(coords, psfmode, fEm=fEm, spectralMethod=spectralMethod)

    #actual stack
    stacked_im  = _stack_stack(method, coords)

    #write image output
    _write_stacked_image(outfile, stacked_im,
                         coords.imagenames[0], stampsize, chanwidth, fEm)

    #plt.figure()
    if plotIt:
        import matplotlib.pyplot as plt
        #NB: si unusedFrequencies = len(coords) on a une division par zero !!!
        #MAIS : ca voudrait aussi dire que pour aucune des iamges on avait de la data dans ce freq
        try:
            stackedImIncludingBorderEffects=stacked_im[int(stampsize/2), int(stampsize/2),0,:]*(len(coords)/(len(coords)-unusedFrequencies))
        except ZeroDivisionError:
            print 'ERORR : your chanwidth is too large and no frequencies were found for ANY of your images in one or more bins'



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
        print data.shape
    else:
        data = 0.*data

def _calculate_sigma2_weights(coords, maskradius=0.):

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
            coord.weight = 0.
        else:
            coord.weight = 1./sigma**2

def _calculate_sigma2_weights_spectral(coords, maskradius=0.):

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
                sigmaF[k]=np.std(tmpdata[:,:,:,k][np.nonzero(tmpdata[:,:,:,k])])

        for k in range(tmpdata.shape[3]):
            tmpdata[:,:,:,k]  = (tmpdata[:,:,:,k]*np.double( X**2+Y**2>maskradius**2))
            sigmaF[k]=np.std(tmpdata[np.nonzero(tmpdata)])
        sigma = np.std(tmpdata[np.nonzero(tmpdata)])
        if sigma == 0:
            coord.weight = 0.
        else:
            coord.weight = 1./sigma**2

def _calculate_amp_weights(coords):
    for i,coord in enumerate(coords):
        amp=data[i,coord.x, coord.y, 0, coord.obsSpecArg]
        coord.weight=1./amp

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
            obsFreq=fEm/(1+coord.z)
            obsFreqArg=int(round((obsFreq-imagesInfo['freq0'][coord.image])/imagesInfo['freqBin'][coord.image]))
        elif spectralMethod=='centralFreq':
            obsFreq=coord.z
            obsFreqArg=int(round((obsFreq-imagesInfo['freq0'][coord.image])/imagesInfo['freqBin'][coord.image]))
        elif spectralMethod=='channel':
            obsFreqArg=coord.z


        #coord.obsSpecArg=obsFreqArg
        if coord.obsSpecArg!=0:
            obsFreqArg=coord.obsSpecArg
        blcx = int(coord.x - stampsize/2 + 0.5)
        blcy = int(coord.y - stampsize/2 + 0.5)
        blcf = obsFreqArg - int(chanwidth/2)
        numberOfZeroLeftF=0

        trcx = blcx + stampsize
        trcy = blcy + stampsize
        trcf = blcf + chanwidth
        numberOfZeroRightF=0
        #coords[i].setWeightToX(chanwidth,X=1)
        if blcf<0:
            listWeights=True
            numberOfZeroLeftF=blcf
            coords[i].setZeroWeightLeft(numberOfZeroLeftF,chanwidth)
            #print 'io'
            blcf=0

        if trcf>imagesInfo['freqlen'][coord.image]-1:
            listWeights=True
            numberOfZeroRightF=trcf-imagesInfo['freqlen'][coord.image]+1
            coords[i].setZeroWeightRight(numberOfZeroRightF, chanwidth)
            trcf=imagesInfo['freqlen'][coord.image]-1

# Currently including positions on the edge of the skymap (in space) will raise an error.
# This could certainly be improved upon.
        if (interval[blcx, trcx] in interval[0, imagesInfo['imagesizes'][coord.image][0]-1]
                and interval[blcy, trcy] in interval[0, imagesInfo['imagesizes'][coord.image][1]-1]):

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
            print blcx, trcx, blcy, trcy, imagesInfo['imagesizes'][coord.image][0]-1
            errorFile=open('errorLog.txt','a')
            errorFile.write(str(coord.image))
            errorFile.close()
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
                coord.weight=([coord.weight for i in range(chanwidth)])

    if Return:
        return data
    #print ([coord.weight for coord in coords])

    #if numberOfZeros is not np.zeros((len(coords), 2)):
    #    for i in range(numberOfZeros.shape[0]):
    #        coords[i].setZeroWeightLeft(numberOfZeros[i][0],chanwidth)
    #        coords[i].setZeroWeightRight(numberOfZeros[i][1],chanwidth)

def _write_stacked_image(imagename, pixels, template_image, stampsize, chanwidth, fEm):
    import os
    import shutil
    from taskinit import ia
#     global stampsize
    if os.access(imagename, os.F_OK): shutil.rmtree(imagename)

    ia.open(template_image)
    beam = ia.restoringbeam()
    cs = ia.coordsys()
    ia.done()

    csnew = cs.copy()
    csnew.setreferencevalue([0.]*2, 'dir')
    csnew.setreferencepixel([int(stampsize/2+0.5)]*2, 'dir')
    if fEm=='no':
        centralFreq=cs.referencevalue()['numeric'][3]
    else:
        centralFreq=fEm
    csnew.setreferencevalue(centralFreq, type='spectral')
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
    pixels = np.zeros(data.shape[1:])

    if method == 'median':
        pixels = np.median(data[0:len(coords),:,:,:,:], 0)

    elif method == 'mean':
        try:
            len(coords[0].weight)
            print 'so many weights'
            pixels=myVeryOwnWeightedAverage(data,coords)
        except TypeError:
            #print data
            #print data.shape
            #np.save('coug', np.array( ([coord.weight for coord in coords]) ))
            pixels = np.average(data[0:len(coords),:,:,:,:], axis=0, weights=([coord.weight for coord in coords]))
            #pixels = np.average(data[0:len(coords),:,:,:,:], axis=0)

    return pixels


def randomizeCoords(coords, beam):
    import random
    import math

    randomcoords = LineStacker.CoordList(coords.imagenames, coords.coord_type,
                             unit=coords.unit)

    for coord in coords:
        dr = random.uniform(5*beam, 10*beam)
        dphi = random.uniform(0, 2*math.pi)
        x = coord.x + dr*math.cos(dphi)
        y = coord.y + dr*math.sin(dphi)
        '''/1? !!!!!
        en dessous: .append coord.z, il faut randomize le z aussi ?
        peut etre peut etre pas ?
        en tout cas j ai mis ca random
        '''
        randomcoords.append(stacker.Coord(x, y, coord.z, coord.weight, coord.image))

    return randomcoords



def noise_estimator(    coords,
                        nrandom=50,
                        imagenames=[],
                        stampsize=32,
                        method='mean',
                        weighting='sigma2',
                        maskradius=None,
                        psfmode='point',
                        fEm=0,
                        chanwidth=100):
    import LineStacker
    import numpy as np
    from taskinit import ia, qa

    print 'estimating noise'
    raise Exception('il faut une condition si Fem==no')
    ia.open(imagenames[0])
    cs=ia.coordsys()
    center=[cs.referencevalue()['numeric'][0], cs.referencevalue()['numeric'][1]]
    Number_Of_Channels=ia.boundingbox()['trc'][3]+1

    #beam = qa.convert(ia.restoringbeam()['major'], 'rad')['value']
    #NB: since there are multiple channels, the value of the
    #restoringbeam is picked from the central frequency
    beam=qa.convert(ia.restoringbeam(int(Number_Of_Channels/2))['major'], 'rad')['value']
    print beam
    mapSizePixels=ia.boundingbox()['trc'][0]-ia.boundingbox()['blc'][0]
    randomPosRangePix=mapSizePixels/2-stampsize/2
    autorizedSize=randomPosRangePix*abs(cs.increment()['numeric'][0])
    randomPosRangeRad=[ [ center[0]-autorizedSize,
                        center[0]+autorizedSize],
                        [ center[1]-autorizedSize,
                        center[1]+autorizedSize]]

    ia.done()

    _allocate_buffers(  imagenames,
                        stampsize,
                        len(coords)*len(imagenames),
                        chanwidth)
    #print 'jambon'
    #print randomPosRangePix, randomPosRangeRad
    print 'buffer allocated'
    dist=([0 for i in range(nrandom)])
    for i in range(nrandom):
        '''
        /!\ z is not random, which means the
        random coordinate will have the same z as the not random one
        (CF: stacker.randomizeCoords)
        '''
        #random_coords = randomizeCoords(coords, beam=beam)
        '''/!\ 2 =^ should be :
        random_coords = stacker.randomizeCoords(coords, beam=beam)
        mais problem a check sur les bornes du random de dr
        '''
        random_coords = LineStacker.randomizeCoords(    coords,
                                                    beam=beam,
                                                    posRangeRad=randomPosRangeRad)
        random_coords = LineStacker.getPixelCoords(random_coords, imagenames)
        _load_stack(random_coords, psfmode, fEm)

        if method == 'mean' and weighting == 'sigma2':
            pass
            '''
            /!\ a faire ! =v
            PUIS MODIF _stack_stack PCQ G REMOVE LE TRUC DES WEIGTHS BB
            '''
            #random_coords = _calculate_sigma2_weights(random_coords, maskradius)
        elif method == 'mean' and weighting == 'sigma':
            pass
            '''
            /!\ a faire ! =v
            '''
            #random_coords = _calculate_sigma_weights(random_coords, maskradius)

        stacked_im  = _stack_stack(method, random_coords)
        dist[i]=(stacked_im[int(stampsize/2+0.5), int(stampsize/2+0.5),0,:])

    return [np.std(dist), np.mean(dist), dist]
