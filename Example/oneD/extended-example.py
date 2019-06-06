import numpy as np
import LineStacker.OneD_Stacker
import LineStacker.analysisTools
import LineStacker.tools.fit as fitTools


#==================================================
#
#                       Stack
#
#==================================================
numberOfSpectra=100
numberOfSpectra2=10

allImages=([0 for i in range(numberOfSpectra+numberOfSpectra2)])
#allImages=([0 for i in range(numberOfSpectra)])
for i in range(numberOfSpectra):
    tempSpectra=np.loadtxt('data/spectra_'+str(i))
    tempCenter=np.loadtxt('data/central_velocity_'+str(i))
    allImages[i]=LineStacker.OneD_Stacker.Image(spectrum=tempSpectra, centralVelocity=tempCenter)


for j in range(numberOfSpectra2):
    tempSpectra=np.loadtxt('data2/spectra_'+str(j))
    tempCenter=np.loadtxt('data2/central_velocity_'+str(j))
    allImages[i+j+1]=LineStacker.OneD_Stacker.Image(spectrum=tempSpectra, centralVelocity=tempCenter)


stacked=LineStacker.OneD_Stacker.Stack(allImages, center='velCenter')

bootstrapped=LineStacker.analysisTools.bootstrapingOneD(allImages, save='amp', center='velCenter', nRandom=10000)
import matplotlib.pyplot as plt
fig=plt.figure()
ax=fig.add_subplot(111)
ax.hist(bootstrapped, bins=20)
#histo=np.histogram(bootstrapped, bins=20)
#fitos=fitTools.GaussFit(fctToFit=histo[0])
#histoBin=histo[1][1]-histo[1][0]
ax.plot(histo[1][1:]-histoBin/2.,fitos)
fig.show()
