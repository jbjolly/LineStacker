"""
**LineStacker.tools.fit Module:**\n
Basic custom fit module for LineStacker.
"""
from scipy.optimize import curve_fit
import numpy as np

c=299792458

def Gauss(  f=[], #this should be a numpy array
            f0=0, #center of your gaussian
            dV='default', #velocity, used to calculate frequency width, set to 'default' to avoid errors
            emittedF=0, #emitted frequency, used to calculate frequency width
            sPeak=1., #amplitude
            dF=1,
            noise=False,
            sigmaNoise=0 ): #frequency width

    if not dF: #if dF is not given by the user it is calculated
        dF=float(dV)*emittedF/c
    sigma=dF/2.3548
    if noise:
        return sPeak*np.exp(-((f-f0)*(f-f0))/(2*sigma*sigma))+np.random.normal(0, sigmaNoise, len(f))
    else:
        return sPeak*np.exp(-((f-f0)*(f-f0))/(2*sigma*sigma))

def gaussFct(x,a,x0,sigma, constant=0):
    """
        Basic Gaussian function
    """
    #+ offset
    if type(x)!=np.ndarray:
        x=np.array(x)

    return constant+a*np.exp(-(x-x0)**2/(2.*sigma**2))
#return float(a)*np.exp(-((x-float(x0))*(x-float(x0)))/(2.*float(sigma)*float(sigma)))

def DoubleGauss(    f=[],
                    amp1=0,
                    amp2=0,
                    f01=0,
                    f02=0,
                    sigma1=0,
                    sigma2=0):
    """
        Sum of two Gaussian function
    """
    return gaussFct(f,amp1,f01,sigma1)+gaussFct(f,amp2,f02,sigma2)

#same amp, same width, centered...
def DoubleGaussPeak(f=[],
                    amp=0,
                    dist=0,
                    sigma=0):

    return gaussFct(f,amp, (f[-1]/2.)-dist/2.,sigma)+gaussFct(f,amp,(f[-1]/2.)+dist/2.,sigma)

def GaussFit(	fctToFit=[],
		fullFreq=[],
		sigma='default',
		returnInfos=False,
        returnError=False):
    """
        Simple Gaussian fitting function. Return order: fitting function, (fitting parameters), (error on fitting parameters)

        Parameters
        ----------
        fctToFit
            One dimensionnal array to be fitted
        fullFreq
            Spectral dimension array. **Not requiered**
        sigma
            initial width for fitting, **'default'** is 3 bins
        returnInfos
            If set to **True** function returns fitted parameters
        returnError
            If set to **True** function returns error on fitted parameters
    """
    if fullFreq==[]:
        fullFreq=range(len(fctToFit))
        if sigma=='default':
            sigma=3
    elif sigma=='default':
        sigma=abs(fullFreq[0]-fullFreq[1])*3
    maxAmp=np.max(fctToFit)
    if type(fctToFit)==np.ndarray:
        maxAmpIndex=fctToFit.tolist().index(maxAmp)
    else:
        maxAmpIndex=fctToFit.index(maxAmp)
    popt, pcov=curve_fit(gaussFct, fullFreq, fctToFit, p0=(maxAmp, fullFreq[maxAmpIndex],sigma))


    if returnInfos:
        if returnError:
            return gaussFct(   fullFreq, popt[0], popt[1], popt[2]), popt, pcov
        else:
            return gaussFct(   fullFreq, popt[0], popt[1], popt[2]), popt
        '''return Gauss(   f=fullFreq,
                    		f0=popt[1],
                    		sPeak=popt[0],dF=popt[2]*2.3548), popt'''
    else:
        if returnError:
            return gaussFct(   fullFreq, popt[0], popt[1], popt[2]), pcov
        else:
            return gaussFct(   fullFreq, popt[0], popt[1], popt[2])


def DoubleGaussFit(         fctToFit=[],
                            fullFreq=[],
                            sigma1='default',
                            sigma2='default',
                            ampScale=0.1,
                            returnAllComp=False,
                            returnInfos=False,
                            returnError=False):
    """
        Double Gaussian (sum of two Gaussians) fitting function. Return order: fitting function, (first Gaussian, second Gaussian), (fitting parameters), (error on fitting parameters)

        Parameters
        ----------
        fctToFit
            One dimensionnal array to be fitted.
        fullFreq
            Spectral dimension array. **Not requiered**
        sigma1
            Initial width of the first Gaussian component for fitting, **'default'** is 3 bins.
        sigma2
            Initial width of the second Gaussian component for fitting, **'default'** is 6 bins.
        ampScale
            Initial amp ratio (between the two Gaussian components) for fitting, **default** is 0.1.
        returnAllComp
            If set to **True** function returns both individual components.
        returnInfos
            If set to **True** function returns fitted parameters.
        returnError
            If set to **True** function returns error on fitted parameters.
    """
    if fullFreq==[]:
        fullFreq=range(len(fctToFit))
        if sigma1=='default':
            sigma1=3
        if sigma2=='default':
            sigma2=6

    else:
        if sigma1=='default':
            sigma1=abs(fullFreq[0]-fullFreq[1])*3
        if sigma2=='default':
            sigma2=abs(fullFreq[0]-fullFreq[1])*6

    maxAmp=np.max(fctToFit)
    if type(fctToFit)==np.ndarray:
    	maxAmpIndex=fctToFit.tolist().index(maxAmp)
    else:
    	maxAmpIndex=fctToFit.index(maxAmp)
    popt,pcov=curve_fit(    DoubleGauss,
                            fullFreq,
                            fctToFit,
                            p0=(    maxAmp,
                                    maxAmp*ampScale,
                                    fullFreq[maxAmpIndex],
                                    fullFreq[maxAmpIndex],
                                    sigma1,
                                    sigma2))
    '''                            p0=(    amp1=maxAmp,
                                    amp2=maxAmp*ampScale,
                                    f01=fullFreq[maxAmpIndex],
                                    f02=fullFreq[maxAmpIndex],
                                    sigma1=sigma1,
                                    sigma2=sigma2))
    '''

    if returnInfos and returnAllComp:
        if returnError:
            return DoubleGauss( f=fullFreq,
                            amp1=popt[0],
                            amp2=popt[1],
                            f01=popt[2],
                            f02=popt[3],
                            sigma1=popt[4],
                            sigma2=popt[5]), gaussFct(fullFreq, popt[0], popt[2], popt[4]), gaussFct(fullFreq, popt[1], popt[3], popt[5]), popt, pcov
        else:
            return DoubleGauss( f=fullFreq,
                            amp1=popt[0],
                            amp2=popt[1],
                            f01=popt[2],
                            f02=popt[3],
                            sigma1=popt[4],
                            sigma2=popt[5]), gaussFct(fullFreq, popt[0], popt[2], popt[4]), gaussFct(fullFreq, popt[1], popt[3], popt[5]), popt
    elif returnInfos:
        if returnError:
            return DoubleGauss( f=fullFreq,
                                amp1=popt[0],
                                amp2=popt[1],
                                f01=popt[2],
                                f02=popt[3],
                                sigma1=popt[4],
                                sigma2=popt[5]), popt, pcov
        else:
            return DoubleGauss( f=fullFreq,
                                amp1=popt[0],
                                amp2=popt[1],
                                f01=popt[2],
                                f02=popt[3],
                                sigma1=popt[4],
                                sigma2=popt[5]), popt
    elif returnAllComp:
        if returnError:
            return DoubleGauss( f=fullFreq,
                                amp1=popt[0],
                                amp2=popt[1],
                                f01=popt[2],
                                f02=popt[3],
                                sigma1=popt[4],
                                sigma2=popt[5]), gaussFct(fullFreq, popt[0], popt[2], popt[4]), gaussFct(fullFreq, popt[1], popt[3], popt[5]), pcov
        else:
            return DoubleGauss( f=fullFreq,
                                amp1=popt[0],
                                amp2=popt[1],
                                f01=popt[2],
                                f02=popt[3],
                                sigma1=popt[4],
                                sigma2=popt[5]), gaussFct(fullFreq, popt[0], popt[2], popt[4]), gaussFct(fullFreq, popt[1], popt[3], popt[5])
    else:
        if returnError:
            return DoubleGauss( f=fullFreq,
                                amp1=popt[0],
                                amp2=popt[1],
                                f01=popt[2],
                                f02=popt[3],
                                sigma1=popt[4],
                                sigma2=popt[5]), pcov
        else:
            return DoubleGauss( f=fullFreq,
                                amp1=popt[0],
                                amp2=popt[1],
                                f01=popt[2],
                                f02=popt[3],
                                sigma1=popt[4],
                                sigma2=popt[5])
def DoublePeakFit(         fctToFit=[],
                            fullFreq=[],
                            returnInfos=False):

    if fullFreq==[]:
        fullFreq=range(len(fctToFit))
        p0Dist=4
        p0Sigma=4
    else:
        p0Dist=4*(fullFreq[-1]-fullFreq[-2])
        p0Sigma=4*(fullFreq[-1]-fullFreq[-2])
    maxAmp=np.max(fctToFit)
    if type(fctToFit)==np.ndarray:
    	maxAmpIndex=fctToFit.tolist().index(maxAmp)
    else:
    	maxAmpIndex=fctToFit.index(maxAmp)
    popt,pcov=curve_fit(    DoubleGaussPeak,
                            fullFreq,
                            fctToFit,
                            p0=(    maxAmp,
                                    p0Dist,
                                    p0Sigma))

    if returnInfos:
        return DoubleGaussPeak( f=fullFreq,
                                amp=popt[0],
                                dist=popt[1],
                            sigma=popt[2]),  popt
    else:
        return DoubleGaussPeak( f=fullFreq,
                                amp=popt[0],
                                dist=popt[1],
                            sigma=popt[2])
