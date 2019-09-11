#This files shows a more extended examplpe use of
#one dimmensional stacking with LineStacker.
#Showcasing the bootstrapping and subsampling tools.
import numpy as np
import LineStacker.OneD_Stacker
import LineStacker.analysisTools
import LineStacker.tools.fit as fitTools
import matplotlib.pyplot as plt


#==================================================
#                       Stack
#==================================================

#In this example all spectra are stored in the dataExtended folder
#There are 50 spectra, 40 having a line amplitude of 1mJy (number 0 to 39)
#while 10 have a line amplitude of 10mJy (number 40 to 49).
#Spectra are simply called spectra_NUMBER with NUMBER going from 0 to 49.
#In addition to each spectra is associated
#an other file called central_velocity_NUMBER
#containing the associated central velocity.
#
#In this example we will perform a median stack of all the spectra.
#In a second time we will perform a bootstrapping analysis,
#indicating the presence of outlayers.
#Finally we will perform a subsampling analysis
#to identify these outlayers.

numberOfSpectra=50
allNames=[]
allImages=([0 for i in range(numberOfSpectra)])
#We start by initializing all spectra as LineStacker.OneD_Stacker.Image
#Which is requiered to do one dimensional stacking
for i in range(numberOfSpectra):
    tempSpectra=np.loadtxt('dataExtended/spectra_'+str(i))
    tempCenter=np.loadtxt('dataExtended/central_velocity_'+str(i))
    allImages[i]=LineStacker.OneD_Stacker.Image(    spectrum=tempSpectra,
                                                    centralVelocity=tempCenter,
                                                    name='spectrum_'+str(i))
    #While the name handling is not requiered it allows easier
    #treatment and understanding of the subsampling.
    allNames.append(allImages[i].name)

#Here we perform the actuall stack
stacked=LineStacker.OneD_Stacker.Stack(allImages, method='median')
#The stack amplitude is stored for later vizualization
stackAmp=np.max(stacked[0])

#==================================================
#              bootstrapping analysis
#==================================================

bootstrapped=LineStacker.analysisTools.bootstraping_OneD(   allImages,
                                                            save='amp',
                                                            nRandom=100000,
                                                            method='median')

#plotting
import matplotlib.pyplot as plt
fig=plt.figure()
ax=fig.add_subplot(111)
counts, bins, patches=ax.hist(bootstrapped, bins=20)
#for vizualization purposes...
ax.set_xticks(bins+(bins[1]-bins[0])/2)
ax.set_xticklabels(([str(bin)[:4] for bin in bins]))
#A red vertical line indicates
#the value of the amplitude of the original stack
ax.axvline(x=stackAmp, color='red')
ax.set_xlabel('stack amplitude')
ax.set_ylabel('count')
ax.set_title('Bootstrapping analysis')
fig.show()

#==================================================
#              subsampling analysis
#==================================================

subSample=LineStacker.analysisTools.subsample_OneD( allImages,
                                                    nRandom=100000,
                                                    method='median')

#plotting
fig=plt.figure()
ax=fig.add_subplot(111)
meanSub=np.mean(([subSample[i] for i in subSample]))
ax.bar(range(len(subSample)), ([subSample[name]-meanSub for name in allNames]))
ax.set_xticks(np.arange(len(allNames))+0.5)
ax.set_xticklabels(allNames, rotation='vertical')
ax.set_ylabel('Spectrum grade - average grade')
ax.set_title('Subsampling analysis')
fig.show()
