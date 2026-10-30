"""
Gaussian fitting of spectra (single, double, and symmetric double peaks).

The fitting functions ignore non-finite samples (NaN channels) and raise
:class:`FitError` (a RuntimeError) when a fit cannot be done or does not
converge: too few valid samples, a constant spectrum, invalid initial
values, or non-finite fitted parameters.
"""
import numpy as np
from scipy.optimize import curve_fit

c=299792458


class FitError(RuntimeError):
    """A fit could not be done or did not converge."""


_FIT_FAILURES = (RuntimeError, ValueError, TypeError, FloatingPointError, ZeroDivisionError)


def _samples(fctToFit, fullFreq, nparams):
    """Finite (x, y) samples to fit, and the full x axis to evaluate on."""
    y = np.asarray(fctToFit, dtype=float).ravel()
    x = np.arange(len(y), dtype=float) if fullFreq is None or len(fullFreq) == 0 else \
        np.asarray(fullFreq, dtype=float).ravel()
    if x.shape != y.shape:
        raise FitError('the spectral axis has {0} values but the data {1}'.format(len(x), len(y)))
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() <= nparams:
        raise FitError('{0} valid samples for {1} parameters'.format(int(ok.sum()), nparams))
    if np.ptp(y[ok]) == 0:
        raise FitError('the data are constant')
    return x[ok], y[ok], x


def _step(x):
    """Typical channel spacing of the valid samples (median spacing, so a
    missing channel does not inflate it)."""
    d = np.diff(np.sort(x))
    d = d[np.isfinite(d) & (d > 0)]
    if d.size == 0:
        raise FitError('the spectral axis has no spacing')
    return float(np.median(d))


def _curve_fit(func, x, y, p0):
    try:
        popt, pcov = curve_fit(func, x, y, p0=p0)
    except _FIT_FAILURES as err:
        raise FitError('fit failed: {0}'.format(err)) from err
    if not np.all(np.isfinite(popt)):
        raise FitError('fit returned non-finite parameters {0}'.format(popt))
    return popt, pcov

def Gauss(  f=None, #this should be a numpy array
            f0=0, #center of your gaussian
            dV='default', #velocity, used to calculate frequency width, set to 'default' to avoid errors
            emittedF=0, #emitted frequency, used to calculate frequency width
            sPeak=1., #amplitude
            dF=1,
            noise=False,
            sigmaNoise=0 ): #frequency width
    """
        Gaussian line on the axis f: peak sPeak at f0 and FWHM dF; with dF=0
        the FWHM is that of the velocity width dV (m/s) at the frequency
        emittedF, dF = dV emittedF / c. noise=True adds Gaussian noise of
        standard deviation sigmaNoise.
    """
    f = np.asarray([] if f is None else f, dtype=float)
    if not dF: #if dF is not given by the user it is calculated
        dF=float(dV)*emittedF/c
    sigma=dF/2.3548
    if noise:
        return sPeak*np.exp(-((f-f0)*(f-f0))/(2*sigma*sigma))+np.random.normal(0, sigmaNoise, len(f))
    else:
        return sPeak*np.exp(-((f-f0)*(f-f0))/(2*sigma*sigma))

def gaussFct(x,a,x0,sigma, constant=0):
    """
        Gaussian a exp(-(x-x0)**2 / (2 sigma**2)) + constant.
    """
    #+ offset
    if type(x)!=np.ndarray:
        x=np.array(x)

    return constant+a*np.exp(-(x-x0)**2/(2.*sigma**2))
#return float(a)*np.exp(-((x-float(x0))*(x-float(x0)))/(2.*float(sigma)*float(sigma)))

def DoubleGauss(    f=None,
                    amp1=0,
                    amp2=0,
                    f01=0,
                    f02=0,
                    sigma1=0,
                    sigma2=0):
    """
        Sum of two Gaussians (amplitudes amp1, amp2, centres f01, f02, widths sigma1, sigma2).
    """
    return gaussFct(f,amp1,f01,sigma1)+gaussFct(f,amp2,f02,sigma2)

#same amp, same width, centered...
def DoubleGaussPeak(f=None,
                    amp=0,
                    dist=0,
                    sigma=0):
    """
        Two Gaussians of the same amplitude and width, dist apart, centred on
        the middle of the axis f.
    """
    f = np.asarray(f, dtype=float)
    mid = 0.5 * (f[0] + f[-1])
    return gaussFct(f,amp, mid-dist/2.,sigma)+gaussFct(f,amp,mid+dist/2.,sigma)

def GaussFit(   fctToFit=None,
                fullFreq=None,
                sigma='default',
                returnInfos=False,
        returnError=False):
    """
        Simple Gaussian fitting function. Return order: fitting function, (fitting parameters), (error on fitting parameters).
        Non-finite samples are ignored; raises FitError (a RuntimeError) if the fit fails.

        Parameters
        ----------
        fctToFit
            One dimensional array to be fitted
        fullFreq
            Spectral dimension array. **Not required**
        sigma
            initial width for fitting, **'default'** is 3 bins
        returnInfos
            If set to **True** function returns fitted parameters
        returnError
            If set to **True** function returns error on fitted parameters
    """
    x, y, fullFreq = _samples(fctToFit, fullFreq, 3)
    if sigma=='default':
        sigma=_step(x)*3
    k = int(np.argmax(y))
    popt, pcov = _curve_fit(gaussFct, x, y, (y[k], x[k], sigma))


    if returnInfos:
        if returnError:
            return gaussFct(   fullFreq, popt[0], popt[1], popt[2]), popt, pcov
        else:
            return gaussFct(   fullFreq, popt[0], popt[1], popt[2]), popt
        # Alternative: return Gauss(f=fullFreq, f0=popt[1],
        #                           sPeak=popt[0], dF=popt[2]*2.3548), popt
    else:
        if returnError:
            return gaussFct(   fullFreq, popt[0], popt[1], popt[2]), pcov
        else:
            return gaussFct(   fullFreq, popt[0], popt[1], popt[2])


def DoubleGaussFit(         fctToFit=None,
                            fullFreq=None,
                            sigma1='default',
                            sigma2='default',
                            ampScale=0.1,
                            returnAllComp=False,
                            returnInfos=False,
                            returnError=False):
    """
        Double Gaussian (sum of two Gaussians) fitting function. Return order: fitting function, (first Gaussian, second Gaussian), (fitting parameters), (error on fitting parameters).
        Non-finite samples are ignored; raises FitError (a RuntimeError) if the fit fails.

        Parameters
        ----------
        fctToFit
            One dimensional array to be fitted.
        fullFreq
            Spectral dimension array. **Not required**
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
    x, y, fullFreq = _samples(fctToFit, fullFreq, 6)
    if sigma1=='default':
        sigma1=_step(x)*3
    if sigma2=='default':
        sigma2=_step(x)*6
    k = int(np.argmax(y))
    popt, pcov = _curve_fit(DoubleGauss, x, y, (y[k], y[k]*ampScale, x[k], x[k], sigma1, sigma2))

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
                                
def DoublePeakFit(         fctToFit=None,
                            fullFreq=None,
                            returnInfos=False):
    """
        Fit of two Gaussians of the same amplitude and width placed
        symmetrically about the middle of the spectral axis (see
        DoubleGaussPeak). Returns the model (and the parameters amp, dist,
        sigma with returnInfos=True). Raises FitError if the fit fails.
    """
    x, y, fullFreq = _samples(fctToFit, fullFreq, 3)
    p0Dist = p0Sigma = 4 * _step(x)
    if not np.all(np.isfinite(fullFreq[[0, -1]])):
        raise FitError('the first and last values of the spectral axis must be finite')
    mid = 0.5 * (fullFreq[0] + fullFreq[-1])

    def model(f, amp, dist, sigma):          # centred on the middle of the full axis
        return gaussFct(f, amp, mid - dist / 2., sigma) + gaussFct(f, amp, mid + dist / 2., sigma)
    popt, pcov = _curve_fit(model, x, y, (np.max(y), p0Dist, p0Sigma))

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
