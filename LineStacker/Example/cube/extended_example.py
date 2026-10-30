import sys
import glob

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

import LineStacker
import LineStacker.image
import LineStacker.analysisTools
import LineStacker.tools.fit as fitTools


def main():
    #coordinates files are selected using the GUI,
    #or all *_coord.txt files of a folder given on the command line
    if len(sys.argv) > 1:
        coordNames=sorted(glob.glob(sys.argv[1].rstrip('/')+'/*_coord.txt'))
    else:
        coordNames=LineStacker.readCoordsNamesGUI()
    coords=LineStacker.readCoords(coordNames)

    #image names are identical to coordinates files,
    #with '.image' replacing '_coord.txt'
    imagenames=([coord.removesuffix('_coord.txt')+'.image' for coord in coordNames])

    #because redshift is used to identify the line center,
    #the emission frequency is also provided.
    #A name not ending in .fits writes a CASA image.
    stacked=LineStacker.image.stack(   coords,
                                            imagenames=imagenames,
                                            fEm=1897420620253.1646,
                                            N_chans=16,
                                            stampsize=16,
                                            outfile='stackResult.image')

    #showing every spectra
    #(spectra are extracted from the central 5x5 pixels of each cube)
    fig=plt.figure(figsize=(10, 16))
    for (i,image_name) in enumerate(imagenames):
        #readCube reads FITS files and CASA images as (x, y, stokes, channel) arrays
        pix=LineStacker.image.readCube(image_name)
        xlen=pix.shape[0]
        spectrum=np.sum(pix[int(xlen/2)-2:int(xlen/2)+3,
                            int(xlen/2)-2:int(xlen/2)+3,
                            :,
                            :], axis=(0,1,2))
        ax=fig.add_subplot(10,5,i+1)
        ax.step(range(len(spectrum)), spectrum, where='mid')
    fig.show()

    #the stack results are stored in the cube stackResult.image
    #Here we open and extract the stacked spectra from the stacked cube,
    #and fit it with 1 and 2 Gaussians
    #Both the spectrum and the fits are then plotted.
    stackResultIm=LineStacker.image.readCube('stackResult.image')
    integratedImage=np.sum(stackResultIm, axis=(0,1,2))

    fited=fitTools.GaussFit(fctToFit=integratedImage)
    doubleFited=fitTools.DoubleGaussFit(fctToFit=integratedImage, returnAllComp=True)

    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.step(range(len(integratedImage)), integratedImage,  'k', label='stack', where='mid')
    ax.plot(fited, 'g', label='Single Gaussian fit')
    ax.plot(doubleFited[0],'r--', label='Double Gaussian Fit, all')
    ax.plot(doubleFited[1],'m', label='Double Gaussian Fit, 1st component')
    ax.plot(doubleFited[2],'c', label='Double Gaussian Fit, 2nd component')
    ax.legend()
    fig.show()

    #Cubes are rebinned,
    #here the linewidths used for rebenning
    #are automatically identified using Gaussian fitting
    rebinnedImageNames=LineStacker.analysisTools.rebin_CubesSpectra(    coords,
                                                                    imagenames)

    #Rebinned cubes are then stacked (written as FITS this time)
    rebinnedStack=LineStacker.image.stack( coords,
                                                imagenames=rebinnedImageNames,
                                                fEm=1897420620253.1646,
                                                N_chans=16,
                                                stampsize=16,
                                                outfile='stackResult_rebinned.fits')
    #Similarly to the previous stack, spectrum is extracted, fited and plotted.
    stackResultIm=LineStacker.image.readCube('stackResult_rebinned.fits')
    integratedImage=np.sum(stackResultIm, axis=(0,1,2))

    #the spectral axis of a stack is in the rest frame: CRVAL3 is the rest
    #frequency of the line, at channel CRPIX3; radio velocities in km/s
    header=fits.getheader('stackResult_rebinned.fits')
    c=299792.458#in km/s
    chans=np.arange(len(integratedImage))-(header['CRPIX3']-1)
    allVel=-c*chans*header['CDELT3']/header['CRVAL3']

    fited=fitTools.GaussFit(fctToFit=integratedImage)
    doubleFited=fitTools.DoubleGaussFit(fctToFit=integratedImage, returnAllComp=True)

    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.step(allVel, integratedImage, 'k', label='stack', where='mid')
    ax.plot(allVel, fited, 'g', label='Single Gaussian fit')
    ax.plot(allVel, doubleFited[0],'r--', label='Double Gaussian Fit, all')
    ax.plot(allVel, doubleFited[1],'m', label='Double Gaussian Fit, 1st component')
    ax.plot(allVel, doubleFited[2],'c', label='Double Gaussian Fit, 2nd component')
    ax.set_xlabel('velocity (km/s)')
    ax.set_title('rebin')
    ax.legend()
    fig.show()

    #Bootstrapping the stack (all CPU cores)
    boot=LineStacker.analysisTools.bootStraping_Cube(   coords,
                                                        imagenames=imagenames,
                                                        fEm=1897420620253.1646,
                                                        N_chans=16,
                                                        stampsize=16,
                                                        nRandom=1000,
                                                        save='ampAndWidth',
                                                        n_jobs=-1)
    return stacked, rebinnedStack, boot


if __name__ == '__main__':
    main()
    plt.show()
