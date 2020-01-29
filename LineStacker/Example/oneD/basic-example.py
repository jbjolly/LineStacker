#This files shows basic examplpe use of
#one dimmensional stacking with LineStacker
import numpy as np
import LineStacker.OneD_Stacker

#In this example we stack 100 spectra that are located in the data folder
#and are named spectra_'+str(i) for i in range(100),
#the lines are identified with a central velocity,
#which can be found in the file 'data/central_velocity_'
#lines can be idenfified in many ways however,
#see LineStacker.OneD_Stacker.Stack for more information
numberOfSpectra=100
allImages=([0 for i in range(numberOfSpectra)])
for i in range(numberOfSpectra):
    tempSpectra=np.loadtxt('data/spectra_'+str(i))
    tempCenter=np.loadtxt('data/central_velocity_'+str(i))
    #initializing all spectra as Image class,
    #this is necessary to use OneD_Stacker.Stack
    allImages[i]=LineStacker.OneD_Stacker.Image(spectrum=tempSpectra, centralVelocity=tempCenter)
stacked=LineStacker.OneD_Stacker.Stack(allImages)
#plotting
import matplotlib.pyplot as plt
fig=plt.figure()
ax=fig.add_subplot(111)
ax.plot(stacked[0])
ax.set_xlabel('Bin number')
ax.set_ylabel('Ampitude')
ax.set_title('Stack of '+str(numberOfSpectra)+' spectra')
fig.show()
