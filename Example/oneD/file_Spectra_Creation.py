import numpy as np
import LineStacker.tools.fit as fitTools

numberOfSpectra=100
lineAmplitude=1
noiseAmplitude=1
lineWidth=200
minVel=-1000
maxVel=1000
numberOfChannels=200
allVelocities=np.linspace(minVel, maxVel, numberOfChannels)

allImages=([0 for i in range(numberOfSpectra)])
for i in range(numberOfSpectra):
    thisLineCenter=np.random.randint(-lineWidth, lineWidth)
    thisLineAmp=abs(np.random.normal(lineAmplitude, lineAmplitude/10.))
    thisLineWidth=abs(np.random.normal(lineWidth, lineWidth/5))
    tempSpectra=fitTools.gaussFct(allVelocities, thisLineAmp, thisLineCenter, thisLineWidth)
    tempSpectra+=np.random.normal(0,noiseAmplitude, len(allVelocities))
    allImages[i]=LineStacker.OneD_Stacker.Image(spectrum=np.array([allVelocities, tempSpectra ]), centralVelocity=thisLineCenter)
    np.savetxt('data/spectra_'+str(i),allImages[i].spectrum)
    centralFile=open('data/central_velocity_'+str(i), 'w')
    centralFile.write(str(allImages[i].centralVelocity))
    centralFile.close()
