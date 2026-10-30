"""Previous (loop-based) OneD Stack, kept as a reference for the tests.
Copied from LineStacker 1.x; only np.NaN -> np.nan changed."""
import numpy as np

def Stack(  Images,
            chansStack='full',
            method='mean',
            center='lineCenterIndex',
            regridFromZ=False,
            regridMethod='scaleToMin'):

    """
        Main (one dimmensional) stacking function.

        Requieres list of Image objects.

        Parameters
        ---------
        Images
            List of images, images have to objects of the Image class (LineStacker.OneD_Stacker.Image).
        chansStack
            Number of channels to stack. Set to 'full' to stack all channels from all images. User input (int) otherwise
        method
            stacking method, **'mean'** and **'median'** supported
        center
            Method to find central frequency of the stack, possible values are:\n
            **"center"**, to stack all spectra center to center,\n
            **'fit'** to use gaussian fitting on the spectrum to determine line center,\n
            **'zero_vel'** to stack on velocity=0 bin,\n
            **'lineCenterIndex'** use the line center initiated with the image,\n
            Or dirrectly defined by the user (int)
        regridFromZ
            if set to True spectra will be regridded to take into account the redshift difference of the different sources. See LineStacker.analysisTools.regridFromZ1D for a more complete description.
        regridMethod
            Used if regridFromZ is True. Can be set either to **'scaleToMin'** or **'scaleToMax'**.
            In the first case all spetcra are regrided to match the smallest redshift (over-gridding), all spectra are regridded to match the highest redshift in the other case (under-gridding).

    """

    if regridFromZ:
        import LineStacker.analysisTools
        Images=LineStacker.analysisTools.regridFromZ1D(Images, regridMethod=regridMethod)

    for (i,image) in enumerate(Images):
        #if image.centerIndex is None:
        if center=='fit':
            #if image.fit==True:
            #    centerIndex=image.centerIndex
            #else: #if the image was not fitted when initiated
            import LineStacker.tools.fit
            if len(image.velocities) > 0:
                tempGaussfit=LineStacker.tools.fit.GaussFit(     fctToFit=image.amp,
                                                    fullFreq=image.velocities,
                                                    returnInfos=True)[1]
                tempVelCenter=tempGaussfit[1] #fit to determine center index
                image.centerIndex=int(round((tempVelCenter-image.velocities[0])/image.velBin))
            else:
                tempGaussfit=LineStacker.tools.fit.GaussFit(     fctToFit=image.amp,
                                                    returnInfos=True)[1]
                image.centerIndex=tempGaussfit[1]

        elif center=='zero_vel': #center on v=0km/s
            try:
                if type(image.velocities)!=list:
                    image.centerIndex=image.velocities.tolist().index(0.)
                else:
                    image.centerIndex=image.velocities.index(0.)
            except ValueError: #if the value 0 is not found, search closest
                centerIndex=-image.velocities[0]/image.velBin
                if centerIndex<0 or centerIndex>len(image.velocities):
                    raise Exception('zero velocity not found')
                image.centerIndex=centerIndex

        elif center=='center': #center on the center of the array
            centerIndex=int(len(image.amp)/2)
            image.centerIndex=centerIndex


        elif center=='lineCenterIndex':
            if image.centerIndex is None:
                raise Exception('center set to lineCenterIndex but no center index is found for line '+str(i))
            else:
                pass
        else: #user defined center
            image.centerIndex=center
            if regridFromZ:
                if regridMethod=='scaleToMin':
                    zToScale=min([tempIm.z for tempIm in Images])
                elif regridMethod=='scaleToMax':
                    zToScale=max([tempIm.z for tempIm in Images])
                zRatio=(1.+image.z)/(1.+zToScale)
                image.centerIndex=int(center*zRatio)

    if chansStack!='full' and type(chansStack)!=int:
        raise Exception("chansStack must be an int or 'full' for full spectra" )

    if chansStack=='full':
        if len(image.velocities) > 0:
            max_size=[int(max([image.centerIndex for image in Images])),
                int(max([len(image.velocities)-image.centerIndex for image in Images]))] #determine the size of the stack output
        elif len(image.frequencies) > 0:
            max_size=[int(max([image.centerIndex for image in Images])),
                int(max([len(image.frequencies)-image.centerIndex for image in Images]))] #determine the size of the stack output

        toStack=np.zeros(( len(Images), max_size[0]+max_size[1])) #create an array, to be filled that will be stacked
        for index in range(-max_size[0],max_size[1]):
            for (i,image) in enumerate(Images):
                if 0<=image.centerIndex+index<len(image.amp): #add to stack only if stacked spectrum defined at this distance from center
                    toStack[i][index+max_size[0]]=image.amp[image.centerIndex+index]
                else:
                    toStack[i][index+max_size[0]]=np.nan
        stacked=np.zeros(max_size[0]+max_size[1]) #initialize the actual stack output
        numberOfSourcePerStack=np.zeros(max_size[0]+max_size[1]) #initialize the number of sources stacked in every channel

        if method=='mean':
            for j in range(max_size[0]+max_size[1]):
                notEmptyPos=np.where( ([not(np.isnan(image[j])) for image in toStack]) )[0] #only stack sources that are not empty at this distance from stack center
                if len(notEmptyPos) > 0:
                    if type(Images[0].weights)!=list and type(Images[0].weights)!=np.ndarray: #if weigths is just one value per spectrum
                        stacked[j]=np.average( ([ toStack[int(ll)][j] for ll in notEmptyPos]), weights=([Images[int(ll)].weights for ll in notEmptyPos]) )
                    else: #if weigth is defined as an array/list (individually for each sectral channel)
                        stacked[j]=np.average( ([ toStack[int(ll)][j] for ll in notEmptyPos]), weights=([Images[int(ll)].weights[j-max_size[0]+Images[ll].centerIndex] for ll in notEmptyPos]) )
                    numberOfSourcePerStack[j]=np.sum([1 for pos in notEmptyPos])

        if method=='median':
            for j in range(max_size[0]+max_size[1]):
                notEmptyPos=np.where( ([not(np.isnan(image[j])) for image in toStack]) )[0]
                if len(notEmptyPos) > 0:
                    stacked[j]=np.median( ([ toStack[int(ll)][j] for ll in notEmptyPos]))
                    numberOfSourcePerStack[j]=np.sum([1 for pos in notEmptyPos])
        return stacked, numberOfSourcePerStack

    elif type(chansStack)==int:
        toStack=np.zeros((len(Images), chansStack))
        for index in range(chansStack):
            for (i, image) in enumerate(Images):
                if len(image.amp)>image.centerIndex-int(chansStack/2.)+index>=0:
                    toStack[i][index]=image.amp[image.centerIndex-int(chansStack/2.)+index]
                else:
                    toStack[i][index]=np.nan
        stacked=np.zeros(chansStack)
        numberOfSourcePerStack=np.zeros(chansStack)

        if method=='mean':
            for j in range(chansStack):
                notEmptyPos=np.where( ([not(np.isnan(image[j])) for image in toStack]) )[0]
                if len(notEmptyPos) > 0:
                    if type(Images[0].weights)!=list and type(Images[0].weights)!=np.ndarray:
                        stacked[j]=np.average( ([ toStack[int(ll)][j] for ll in notEmptyPos]), weights=([Images[int(ll)].weights for ll in notEmptyPos]) )
                    else:
                        stacked[j]=np.average( ([ toStack[int(ll)][j] for ll in notEmptyPos]), weights=([Images[int(ll)].weights[image.centerIndex-int(chansStack/2.+0.5)+j] for ll in notEmptyPos]) )
                    numberOfSourcePerStack[j]=np.sum([1 for pos in notEmptyPos])
        elif method=='median':
            for j in range(chansStack):
                notEmptyPos=np.where( ([not(np.isnan(image[j])) for image in toStack]) )[0]
                if len(notEmptyPos) > 0:
                    stacked[j]=np.median( ([ toStack[int(ll)][j] for ll in notEmptyPos]))
                    numberOfSourcePerStack[j]=np.sum([1 for pos in notEmptyPos])
        return stacked, numberOfSourcePerStack
