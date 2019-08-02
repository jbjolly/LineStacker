import sys
import numpy as np
"""
	Library including statistical tools for stacking
"""

def bootstraping_OneD(	Images,
						nRandom=1000,
						chansStack='full',
						method='mean',
						center='z',
						velOrFreq='vel',
						save='all'):

	"""
        Performs bootstrapping stack of spectra,
		see stacker.OneD_Stacker.stack for further information on stack parametres
        Parameters
        ---------
		Images
			a list of stacker.OneD_Stacker.Images
		nRandom
			number of boostrap itterations\n
			default is 1000
		chansStack
			number of channels to stack, either a fixed number or 'full' for entire spectra\n
			default is 'full'
		method
			stacking method, 'mean' or 'median'\n
			default is 'mean'
		center
			method to center spectra\n
			possible values are 'z', 'fit', 'zero_vel', 'center' or user input (must be int)\n
			see stacker.OneD_Stacker.stack for further information on centering methods
			default is 'z'
		velOrFreq
			velocity or frequency mode (chose according to your spectral axis type)\n
			possible values are 'vel' or 'freq'\n
			default is 'vel'
		save
			data to save at each bootstrap itterations\n
			possible values are 'all', 'amp' and 'ampAndWidth'\n
			'all' saves the full stack at each bootstrap itteration /!\\ caution, can be memory expensive\n
			'amp' saves the amplitude of the line (maximum value of the stack) at each bootstrap itteration, fastest\n
			'ampAndWidth' fits the line with a gaussian and saves the corresponding amplitude and width at each bootstrap itteration, can be cpu exppensive\n
			default is 'all'
    """
	if 'stacker.OneD_Stacker' not in sys.modules:
		import LineStacker.OneD_Stacker

	if not isinstance(Images[0], LineStacker.OneD_Stacker.Image):
		raise Exception('images must belong to the Image class of stacker.OneD_Stacker please initialize your images accordingly')

	baseStack=LineStacker.OneD_Stacker.Stack(	Images,
											chansStack=chansStack,
											method=method,
											center=center,
											velOrFreq=velOrFreq)
	saved=([0 for i in range(nRandom+1)])
	'''/!\
	'''
	import LineStacker.tools.fit


	if save=='all':
		saved[0]=baseStack[0]
	elif save=='amp':
		saved[0]=np.max(baseStack[0])
	elif save=='ampAndWidth':
		if 'LineStacker.tools.fit' not in sys.modules:
			print 'fit module not imported'
			import LineStacker.tools.fit
		fited=LineStacker.tools.GaussFit(fctToFit=baseStack[0], returnInfos=True)
		saved[0]=[fited[1][0],fited[1][2]]

	g=0

	print '\nBootstrap starting ...\n'
	for n in range(nRandom):

		LineStacker.tools.ProgressBar(n, nRandom)

		newImages=([Images[np.random.randint(len(Images))] for i in Images])
		tempStack=LineStacker.OneD_Stacker.Stack(	newImages,
												chansStack=chansStack,
												method=method,
												center=center,
												velOrFreq=velOrFreq)
		if save=='all':
			saved[n+1]=tempStack[0]
		elif save=='amp':
			saved[n+1]=np.max(tempStack[0])
		elif save=='ampAndWidth':
			fited=LineStacker.tools.fit.GaussFit(fctToFit=tempStack[0], returnInfos=True)
			saved[n+1]=[fited[1][0],fited[1][2]]

	return saved

def bootStraping_Cube(	coords,
				        outfile,
				        stampsize = 32,
				        imagenames= [],
				        method = 'mean',
				        weighting = None,
				        maxmaskradius=0,
				        fEm = 0,
				        chanwidth=30,
						nRandom=1000,
						save='amp'):#'all', 'amp', 'ampAndWidth'

	"""
        Performs bootstrapping stack of cubes,
		see stacker.line_image.stack for further information on stack parametres
        Parameters
        ---------
		coords
			A stacker.coordList object of all target coordinates
		outfile
			Target name for stacked image
		stampsize
			size of target image in pixels\n
			default is 32
		imagenames
			Name of images to extract cubes from
		method
			stacking method, 'mean' or 'median'\n
			default is 'mean'
		weighting
			weighting method to use if stacking method is mean\n
			possible values are 'sigma2', 'sigma2F', '1/A', 'None', 1 or user input (float)\n
			see stacker.line-image for a complete description of weighting methods\n
			default is None
		maskradius
			radius of the mask used to blank the centre pixels in weight calculation\n
			default is 0
		fEm
			rest emission frequency of the line
		chanwidth
			number of channels of the resulting stack
		nRandom
			number of boostrap itterations\n
			default is 1000
		save
			data to save at each bootstrap itterations\n
			possible values are 'all', 'amp' and 'ampAndWidth'\n
			'all' saves the full stack at each bootstrap itteration /!\\ caution, can be memory expensive\n
			'amp' saves the amplitude (maximum value of the stack) of the line, at each bootstrap itteration, fastest\n
			'ampAndWidth' fits the line with a gaussian and saves the corresponding amplitude and width at each bootstrap itteration, can be cpu exppensive\n
			'ouflow' fits the line with two gaussian components and saves the stack parameters at each bootstrap itteration, can be cpu exppensive\n
			for 'amp', 'ampAndWidth' and 'ouflow' the line is obtained by summing all pixels inside the stack stamp\n
			default is 'all'

    """

	import LineStacker.line_image
	import LineStacker
	import LineStacker.tools
	if coords.coord_type == 'physical':
		coords = LineStacker.getPixelCoords(coords, imagenames)

	LineStacker.line_image._allocate_buffers(coords.imagenames, stampsize, len(coords), chanwidth)

	if method=='mean' and coords[0].weight==1 and weighting=='sigma2':
		LineStacker.line_image._calculate_sigma2_weights(coords, maskradius=maskradius)
	if method=='mean' and coords[0].weight==1 and weighting=='sigma2F':
		LineStacker.line_image._calculate_sigma2_weights_spectral(coords, maskradius=maskradius)
	if method=='mean' and coords[0].weight==1 and weighting=='1/A':
		LineStacker.line_image._calculate_amp_weights(coords)

    #fills data with skymap accordingly
	data=LineStacker.line_image._load_stack(coords, psfmode, fEm=fEm, Return=True)
	stacked_im  = LineStacker.line_image._stack_stack(method, coords)

	saved=([0 for i in range(nRandom+1)])
	if save=='all':
		saved[0]=stacked_im
	elif save=='amp':
		saved[0]=np.max(np.sum(stacked_im, axis=(1,2,3)))
	elif save=='ampAndWidth':
		import LineStacker.tools.fit
		fited=LineStacker.tools.fit.GaussFit(fctToFit=np.sum(stacked_im, axis=(1,2,3)),returnInfos=True)
		saved[0]=[fited[1][0],fited[1][2]]
	elif save=='outflow':
		import LineStacker.tools.fit
		fited=LineStacker.tools.fit.DoubleGaussFit(fctToFit=np.sum(stacked_im, axis=(1,2,3)),returnInfos=True)
		saved[0]=[fited[1][1],fited[1][5]]
	for n in range(nRandom):
		LineStacker.tools.ProgressBar(n, nRandom)
		newList=np.random.randint(0,len(coords),len(coords))
		newCoords=([coords[i] for i in newList])
		if method=='mean':
			try:
				len(newCoords[0].weight)
				#print 'so many weights'
				pixels=LineStacker.line_image.myVeryOwnWeightedAverage(np.array([data[i] for i in newList]), newCoords)
			except TypeError:
				pixels=np.average([data[i] for i in newList], axis=0, weights=([coord.weight for coord in newCoords]))
		elif method=='median':
			pixels=np.median([data[i] for i in newList], axis=0)

		if save=='all':
			saved[n+1]=pixels
		elif save=='amp':
			saved[n+1]=np.max(np.sum(pixels, axis=(1,2,3)))
		elif save=='ampAndWidth':
			fited=LineStacker.tools.fit.GaussFit(fctToFit=np.sum(pixels, axis=(1,2,3)),returnInfos=True)
			saved[n+1]=[fited[1][0],fited[1][2]]
		elif save=='outflow':
			try:
				fited=LineStacker.tools.fit.DoubleGaussFit(fctToFit=np.sum(pixels, axis=(1,2,3)),returnInfos=True)
				saved[n+1]=[fited[1]]
			except RuntimeError:
				saved[n+1]=[-100 for i in range(6)]

	return saved

def congrid(a, newdims, method='linear', centre=False, minusone=False):
	import scipy.interpolate
	import scipy.ndimage

	if not a.dtype in [np.float64, np.float32]:
	    a = np.cast[float](a)

	m1 = np.cast[int](minusone)
	ofs = np.cast[int](centre) * 0.5
	old = np.array( a.shape )
	ndims = len( a.shape )
	if len( newdims ) != ndims:
	    print "[congrid] dimensions error. " \
	          "This routine currently only support " \
	          "rebinning to the same number of dimensions."
	    return None
	newdims = np.asarray( newdims, dtype=float )
	dimlist = []

	if method == 'neighbour':
	    for i in range( ndims ):
	        base = np.indices(newdims)[i]
	        dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
	                        * (base + ofs) - ofs )
	    cd = np.array( dimlist ).round().astype(int)
	    newa = a[list( cd )]
	    return newa

	elif method in ['nearest','linear']:
	    # calculate new dims
	    for i in range( ndims ):
	        base = np.arange( newdims[i] )
	        dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
	                        * (base + ofs) - ofs )
	    # specify old dims
	    olddims = [np.arange(i, dtype = np.float) for i in list( a.shape )]

	    # first interpolation - for ndims = any
	    mint = scipy.interpolate.interp1d( olddims[-1], a, kind=method )
	    newa = mint( dimlist[-1] )

	    trorder = [ndims - 1] + range( ndims - 1 )
	    for i in range( ndims - 2, -1, -1 ):
	        newa = newa.transpose( trorder )

	        mint = scipy.interpolate.interp1d( olddims[i], newa, kind=method )
	        newa = mint( dimlist[i] )

	    if ndims > 1:
	        # need one more transpose to return to original dimensions
	        newa = newa.transpose( trorder )

	    return newa
	elif method in ['spline']:
	    oslices = [ slice(0,j) for j in old ]
	    oldcoords = np.ogrid[oslices]
	    nslices = [ slice(0,j) for j in list(newdims) ]
	    newcoords = np.mgrid[nslices]

	    newcoords_dims = range(np.rank(newcoords))
	    #make first index last
	    newcoords_dims.append(newcoords_dims.pop(0))
	    newcoords_tr = newcoords.transpose(newcoords_dims)
	    # makes a view that affects newcoords

	    newcoords_tr += ofs

	    deltas = (np.asarray(old) - m1) / (newdims - m1)
	    newcoords_tr *= deltas

	    newcoords_tr -= ofs

	    newa = scipy.ndimage.map_coordinates(a, newcoords)
	    return newa
	else:
	    print "Congrid error: Unrecognized interpolation type.\n", \
	          "Currently only \'neighbour\', \'nearest\',\'linear\',", \
	          "and \'spline\' are supported."
	    return None

def rebin_OneD(images):
	"""
        Rebin a list of images so that all width have the same width as the smallest width\r
		/!\\ Lines must be visible before stacking to operate rebinning
        Parameters
        ---------
		images
			A list of stacker.OneD_Stacker.images
	"""
	import scipy.interpolate
	import scipy.ndimage
	minWidth=1e9
	minWidthIndex=1e9
	allWidth=([0 for i in images])
	if images[0].fit:
		for (i, image) in enumerate(images):
			tempWidth=image.gaussfit[1][2]
			if tempWidth<minWidth:
				minWidth=tempWidth
				minWidthIndex=i
			allWidth[i]=tempWidth
	else:
		import LineStacker.tools.fit
		import matplotlib.pyplot as plt
		fig=plt.figure()
		for (i, image) in enumerate(images):
			tempWidth=LineStacker.tools.fit.GaussFit(fctToFit=image.amp, returnInfos=True)[1][2]
			fitos=LineStacker.tools.fit.GaussFit(fctToFit=image.amp)
			ax=fig.add_subplot(6,5,i+1)
			ax.plot(image.amp)
			ax.plot(fitos,'r')
			if tempWidth<minWidth:
				minWidth=tempWidth
				minWidthIndex=i
			allWidth[i]=tempWidth
		fig.show()
	import copy
	allNewImages=[]

	for (i, image) in enumerate(images):
		newImage=copy.deepcopy(image)
		if i!=minWidthIndex:
			lenRatio=minWidth/allWidth[i]
			for attr in vars(newImage):
				try:
					newLen=int(len(vars(newImage)[attr])*lenRatio)
					if type(vars(newImage)[attr])!=str:
						if attr=='spectrum':
							if vars(newImage)[attr].shape[1]==2:
								vars(newImage)[attr]=congrid(np.array(vars(newImage)[attr]), [newLen, 2])
							elif vars(newImage)[attr].shape[0]==2:
								vars(newImage)[attr]=congrid(np.array(vars(newImage)[attr]), [2, newLen])
						elif attr=='gaussfit':
							pass
						else:
							vars(newImage)[attr]=congrid(np.array(vars(newImage)[attr]), [newLen])
				except TypeError:
					pass
				if attr=='velBin':
					vars(newImage)[attr]=vars(newImage)[attr]/lenRatio
				if attr=='freqbin':
					vars(newImage)[attr]=vars(newImage)[attr]/lenRatio

			if newImage.fit:
				newImage.fitAgain()
		allNewImages.append(newImage)
	print '\nrebinning done\n'
	return allNewImages


def rebin_CubesSpectra(	coords,
						imagenames,
						regionSize=False,
						widths=False):

	"""
        Rebin a list of image-cubes so that all width have the same width as the smallest width\r
		/!\\ Lines must be visible before stacking to operate rebinning\r
		/!\\ Only one coord per image is necessary for rebinning

        Parameters
        ---------
		coords
			A coordList object of all target coordinates.
		imagenames
			Name of images to rebin
		regionSize
			size (in pixels) to extract the spectra from\r
			if set to False, spectra will be extracted solely from the coord pixel
		widths
			widths of the lines\r
			if set to 'False' the spectra will be fitted with a gaussian to extract the width
	"""

	print 'rebinning...\r'
	from taskinit import ia
	#NB:ONE LINE PER CUBE, duh
	import LineStacker.tools
	imageAsArray=([0 for i in coords])
	allCoordsSys=([0 for i in coords])
	if not widths:
		import LineStacker.tools.fit as fitTools
		if coords.coord_type == 'physical':
			coords = LineStacker.getPixelCoords(coords, imagenames)
		allWidths=[]
		minWidth=1e30
		import matplotlib.pyplot as plt
		fig=plt.figure()
		for (i,image) in enumerate(coords.imagenames):
			coord=coords[i]
			ia.open(image)
			imageAsArray[i]=ia.getregion()
			allCoordsSys[i]=ia.coordsys()
			ia.done()
			ax=fig.add_subplot(6,5,i+1)
			if regionSize:
				pixels=imageAsArray[i][coord.x-int(regionSize/2):coord.x+int(regionSize/2):, coord.y-int(regionSize/2):coord.y+int(regionSize/2),0,:]
				spectra=np.sum(pixels, axis=(0,1))
			else:
				spectra=imageAsArray[i][coord.x, coord.y,0,:]
			fitted=fitTools.GaussFit(fctToFit=spectra, returnInfos=True)
			ax.plot(spectra)
			ax.plot(fitted[0], 'r')
			ax.text(0.2, 0.8, str(fitted[1][2])[:6], horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
			allWidths.append(fitted[1][2])
			if allWidths[i]<minWidth:
				minWidth=allWidths[i]
	else:
		minWidth=np.min(widths)
		allWidths=widths
	widthsRatios=minWidth/np.array(allWidths)
	for (i,image) in enumerate(imagenames):
		LineStacker.tools.ProgressBar(i,len(coords))
		if allWidths[i]!=minWidth:
			if widths:
				ia.open(image)
				imageAsArray[i]=ia.getregion()
				allCoordsSys[i]=ia.coordsys()
				ia.done()
			else:
				pass
			newShape=list(imageAsArray[i].shape)
			newShape[3]=int(newShape[3]*widthsRatios[i])
			newShape=[newShape[0],newShape[1],newShape[3]]
			newArray=congrid(imageAsArray[i][:,:,0,:], newShape)

			newArray=newArray.reshape(newShape[0], newShape[1],1, newShape[2])
			oldInc=allCoordsSys[i].increment()['numeric'][3]
			allCoordsSys[i].setincrement(oldInc*imageAsArray[i].shape[3]/newShape[2],type='spectral')
			tempOutname=image.split('.')
			if len(tempOutname)==2:
				outname=tempOutname[0]+'__ULTIMATE_SpectralRebinned.image'
			elif tempOutname>2:
				outname=''
				for nameSplit in tempOutname[:-2]:
					outname+=nameSplit+'.'
				outname+=tempOutname[-2]
				outname+='__ULTIMATE_SpectralRebinned.image'
			else:
				outname=tempOutname+'__ULTIMATE_SpectralRebinned.image'
			print outname
			ia.fromarray(	outname,
							pixels=newArray,
							csys=allCoordsSys[i].torecord(),
							overwrite=True)
			ia.done()
		else:
			ia.open(image)
			tempOutname=image.split('.')
			if len(tempOutname)==2:
				outname=tempOutname[0]+'__ULTIMATE_SpectralRebinned.image'
			elif tempOutname>2:
				outname=''
				for nameSplit in tempOutname[:-2]:
					outname+=nameSplit+'.'
				outname+=tempOutname[-2]
				outname+='__ULTIMATE_SpectralRebinned.image'
			else:
				outname=tempOutname+'__ULTIMATE_SpectralRebinned.image'
			print outname
			ia.subimage(outfile=outname, overwrite=True)
			ia.done()
def randomizeSample(	sample,
						minSize=3,
						maxSize=None):
	if maxSize==None:
		newSampleSize=np.random.randint(minSize, len(sample)+1)
	else:
		newSampleSize=np.random.randint(minSize, maxSize+1)
	newSample=np.array(sample)
	np.random.shuffle(newSample)
	newSample=newSample[:newSampleSize]
	return newSample

def maximizeAmp(spectra):
	return np.max(spectra)

def maximizeSNR(spectra):
	import LineStacker.tools.fit as fitTools
	#try:
	fitos=fitTools.GaussFit(fctToFit=spectra, returnInfos=True)
	amp, center, width=fitos[1]
	#e
	toSTD=[]
	leftLimSpectra=int(center-width*2.35)
	rightLimSpectra=int(center+width*2.35)
	if rightLimSpectra<len(spectra)-1:
		toSTD.append(spectra[rightLimSpectra:][:])
	if leftLimSpectra>0:
		toSTD.append(spectra[:leftLimSpectra][:])
	if len(toSTD)==2:
		toSTD=np.concatenate(toSTD)
	elif len(toSTD)==1:
		toSTD=toSTD[0]
	else:
		return 0
	toSTD=np.array(toSTD).flatten()
	return float(amp)/np.std(toSTD)

def maximizeOutflow(spectra):
	import LineStacker.tools.fit as fitTools
	try:
		doubleFitos=fitTools.DoubleGaussFit(fctToFit=spectra, returnInfos=True)[1]
	except RuntimeError:
		return 0
	if doubleFitos[0]>doubleFitos[1]:
		secondCompIndex=1
		area=doubleFitos[1]*doubleFitos[5]**2
	else:
		secondCompIndex=0
		#area=doubleFitos[0]*doubleFitos[4]
		area=doubleFitos[0]*doubleFitos[4]**2
	return area

def subsample_OneD(	images,
					nRandom=10000,
					maxTest=maximizeSNR,
					**kwargs):
	'''
		Randomly ressamples spectra and grades them accordingly
		to a given grade function

		Parameters
		---------

		images
			A list of stacker.OneD_Stacker.images
		nRandom
			Number of itterations of the Monte-Carlo process
		maxTest
			function test to grade the sources\r
			build your own, or use existing: maximizeAmp, maximizeSNR, maximizeOutflow
		Any other argument for LineStacker.OneD_Stacker.Stack
	'''
	import LineStacker.tools
	import LineStacker.OneD_Stacker
	imagesGrades={}
	for (i,image) in enumerate(images):
		if image.name=='':
			image.name=str(i)
		imagesGrades[str(image.name)]=0
	for n in range(nRandom):
		LineStacker.tools.ProgressBar(n,nRandom)
		newImages=randomizeSample(images)
		tempStack=LineStacker.OneD_Stacker.Stack(newImages, **kwargs)
		testResult=maxTest(tempStack[0])
		for image in newImages:
			imagesGrades[str(image.name)]+=testResult
	return imagesGrades

def stack_estimator(    coords,
                        nRandom=100,
                        imagenames=[],
                        stampsize=1,
                        method='mean',
                        chanwidth=30,
						lowerLimit='default',
						upperLimit='default',
						**kwargs):
    """
        Performs stacks at random positions a set number of times.
		Allows to probe the relevance of stack through stacking random positions
		as a Monte Carlo process


        returns: Estimate of stacked flux assuming point source.

        Parameters
        ---------
        coords
            A coordList object of all target coordinates.
	    nRandom
			Number of itterations
        imagenames
            Name of imagenames to stack
        stampsize
			Size of the stamp to stack (because only the central pixel is extracted anyway, this should be kept to default to enhance efficiency)
    	method
			Method for stacking, see LineStacker.line_image.stack
        chanwidth
			Number of channels in the stack
		lowerLimit
			Lower spatial limit (distance from stacking position) to random new random position.\n
			Default is 5 beams
		upperLimit
			Upper spatial limit (distance from stacking position) to random new random position.\n
			Default is 10 beams
	"""
    import LineStacker
    import numpy as np
    from taskinit import ia, qa

    print 'estimating stack relevance'

    ia.open(imagenames[0])
    cs=ia.coordsys()
    center=[cs.referencevalue()['numeric'][0], cs.referencevalue()['numeric'][1]]
    Number_Of_Channels=ia.boundingbox()['trc'][3]+1

    #NB: since there are multiple channels, the value of the
    #restoringbeam is picked from the central frequency
    beam=qa.convert(ia.restoringbeam(int(Number_Of_Channels/2))['major'], 'rad')['value']
    ia.done()

    _allocate_buffers(  imagenames,
                        stampsize,
                        len(coords),
                        chanwidth)

    dist=([0 for i in range(nRandom)])
    for i in range(nRandom):
        '''
        /!\\ z is not random, which means the
        random coordinate will have the same z as the not random one
        (CF: randomizeCoords)
        '''
		random_coords = randomizeCoords(	coords,
											beam=beam,
											lowerLimit=lowerLimit,
											upperLimit=upperLimit)

        random_coords = LineStacker.getPixelCoords(random_coords, imagenames)
        _load_stack(	random_coords,
						psfmode,
						fEm=fEm,
						spectralMethod=spectralMethod)

        stacked_im  = _stack_stack(method, random_coords)
		#from the stack of random positions only the
		#spectrum of the central pixel is kept
        dist[i]=(stacked_im[int(stampsize/2+0.5), int(stampsize/2+0.5),0,:])

    return [np.std(dist), np.mean(dist), dist]

def randomizeCoords(	coords,
						beam=beam,
						lowerLimit='default',
						upperLimit='default'):
    """
        Function to randomize stacking coordinates, new ramdom position uniformally randomized, centered on the original stacking coordinate

        Parameters
        ---------
        coords
            A coordList object of all target coordinates.
	    beam
			beam size, used for default values of random range
		lowerLimit
			Lower spatial limit (distance from stacking position) to random new random position.\n
			Default is 5 beams
		upperLimit
			Upper spatial limit (distance from stacking position) to random new random position.\n
			Default is 10 beams
	"""
    import numpy as np
    import math
	if lowerLimit=='default':
		lowerLimit=beam*5.
	if upperLimit=='default':
		upperLimit=beam*10.
    randomcoords = LineStacker.CoordList(coords.imagenames, coords.coord_type,
                             unit=coords.unit)

    for coord in coords:
        dr = np.random.uniform(lowerLimit, upperLimit)
        dphi = np.random.uniform(0, 2*math.pi)
        x = coord.x + dr*math.cos(dphi)
        y = coord.y + dr*math.sin(dphi)
        randomcoords.append(stacker.Coord(x, y, coord.z, coord.weight, coord.image))

    return randomcoords
