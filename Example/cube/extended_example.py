import LineStacker
import LineStacker.line_image
import LineStacker.analysisTools
import numpy as np
from taskinit import ia
import LineStacker.tools.fit as fitTools
import matplotlib.pyplot as plt


#coordinates files are selected using the GUI
coordNames=LineStacker.readCoordsNamesGUI()
coords=LineStacker.readCoords(coordNames)

#image names are identical to coordinates files,
#with '.image' replacing '_coords.txt'
imagenames=([coord.strip('_coords.txt')+'.image' for coord in coordNames])

#because redshift is used to identify the line center,
#the emission frequency is also provided
stacked=LineStacker.line_image.stack(   coords,
                                        imagenames=imagenames,
                                        fEm=1897420620253.1646,
                                        stampsize=8)

#showing every spectra
#(spectra are extracted from the central 5x5 pixels of each cube)
fig=plt.figure()
for (i,image) in enumerate(imagenames):
    ia.open(image)
    pix=ia.getchunk()
    ia.done()
    xlen=pix.shape[0]
    spectrum=np.sum(pix[int(xlen/2)-2:int(xlen/2)+3,
                        int(xlen/2)-2:int(xlen/2)+3,
                        :,
                        :], axis=(0,1,2))
    ax=fig.add_subplot(10,5,i+1)
    ax.step(range(len(spectrum)), spectrum, where='mid')
fig.show()

#the stack results are stored in the cube stackResult.image by default
#Here we open and extract the stacked spectra from the stacked cube,
#and fit it with 1 and 2 Gaussians
#Both the spectrum and the fits are then plotted.
ia.open('stackResult.image')
stackResultIm=ia.getchunk()
ia.done()
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

#Rebinned cubes are then stacked
rebinnedStack=LineStacker.line_image.stack(   coords,
              		                          imagenames=rebinnedImageNames,
              		                          fEm=1897420620253.1646,
              		                          stampsize=8)
#Similarly to the previous stack, spectrum is extracted, fited and plotted.
ia.open('stackResult.image')
stackResultIm=ia.getchunk()
ia.done()
integratedImage=np.sum(stackResultIm, axis=(0,1,2))

fited=fitTools.GaussFit(fctToFit=integratedImage)
doubleFited=fitTools.DoubleGaussFit(fctToFit=integratedImage, returnAllComp=True)

fig=plt.figure()
ax=fig.add_subplot(111)
ax.step(range(len(integratedImage)), integratedImage, 'k', label='stack', where='mid')
ax.plot(fited, 'g', label='Single Gaussian fit')
ax.plot(doubleFited[0],'r--', label='Double Gaussian Fit, all')
ax.plot(doubleFited[1],'m', label='Double Gaussian Fit, 1st component')
ax.plot(doubleFited[2],'c', label='Double Gaussian Fit, 2nd component')
ax.set_title('rebin')
ax.legend()
fig.show()
