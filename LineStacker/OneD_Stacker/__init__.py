"""
Stacking of 1D spectra: each spectrum is a :class:`Spectrum`, which knows
where its line is; :func:`Stack` aligns the spectra on their lines and
averages them channel by channel. Velocities are in km/s.
"""
import numpy as np

import LineStacker
import LineStacker.tools.fit

#: speed of light in km/s: velocities of spectra are in km/s
c = 299792.458

VELOCITY_CONVENTIONS = ('radio', 'optical', 'relativistic')

#: Largest distance (in channels) of a spectral value from an evenly spaced
#: axis that Spectrum accepts.
MAX_AXIS_DEVIATION = 0.05


class Spectrum():
    """

        One spectrum to stack: its fluxes, optionally its spectral values, and
        where its line is.

        Parameters
        ---------
        spectrum
            Spectrum should be of the shape [N,2] or [2,N],
            where N is the number of spectral bins.
            The second dimension being the flux
            and the first the spectral values, which must be evenly spaced
            (increasing or decreasing; within 5% of a channel).
        amp
            Instead of defining spectrum one can define only the amplitude (flux).
        velocities
            In the case where amp is defined
            it is still possible to define the velocities
        frequencies
            In the case where amp is defined
            it is still possible to define the frequencies
        z
            Redshift, one of the possible way to define central frequency.
            If using redshift **fEmLine** (emission frequency of the line) should be defined.
            In addition the spectrum should consist of BOTH amplitudes and spectral information (frequencies or velocities), and not amplitude alone.
            Velocities are then in km/s relative to the rest frequency fEmLine (as e.g. a CASA spectrum with restfreq=fEmLine), in the convention
            **velConvention**: the line is at v = c z/(1+z) (radio), c z (optical) or c((1+z)**2-1)/((1+z)**2+1) (relativistic).
        fEmLine
            Emission frequency of the line, needed if redshift argument is used.
        centerIndex
            channel index of the line center, another possible way to define line center.
        centralFrequency
            (observed) Frequency of the line center, another possible way to define line center.
        centralVelocity
            Observed velocity of the line center, another possible way to define line center.
        weights
            The weighting scheme to use.\n
            Can be set to **'1/A'** or **'sigma2'**.
            If **sigma2** is used std of the entire spectra is used UNLESS, **fit** is set to **True**, in which case the central part around the line is excluded from std calculation (central here means one FWHM on each side of the center of the line).\n
            Alternatively user input can be used, float or list (or array)
        name
            name of the spectrum, can allow easier identification
        fit
            If fit is set to **True**, the spectrum will be fitted with a gaussian to try and identify the line center (as well as amplitude if weight is set to 1/A)
        velOrFreq
            Defining if spectral dimension is velocity or frequency.
        velConvention
            Velocity convention of the velocities (km/s), used with z: **'radio'** (default, as CASA and LineStacker.uv),
            **'optical'** or **'relativistic'**.
    """

    def __init__(   self,
                    spectrum=None,
                    amp=None,
                    velocities=None,
                    frequencies=None,
                    z=None,
                    fEmLine=False,
                    centerIndex=None,
                    centralVelocity=None,
                    centralFrequency=None,
                    weights=1,
                    name='',
                    fit=False,
                    velOrFreq='vel',
                    velConvention='radio'):

        spectrum = np.asarray([] if spectrum is None else spectrum, dtype=float)
        amp = [] if amp is None else amp
        velocities = [] if velocities is None else velocities
        frequencies = [] if frequencies is None else frequencies
        if velOrFreq not in ('vel', 'freq'):
            raise ValueError("velOrFreq must be 'vel' or 'freq'")
        if velConvention not in VELOCITY_CONVENTIONS:
            raise ValueError('velConvention must be one of {0}'.format(VELOCITY_CONVENTIONS))
        self.velConvention = velConvention
        self.velOrFreq=velOrFreq
        self.name=name
        self.z=z
        self.fEmLine=fEmLine
        self.fit=fit
        self.centerIndex=centerIndex
        self.centralVelocity=centralVelocity
        self.centralFrequency=centralFrequency
        self.spectrum=spectrum
        self.velocities=velocities
        self.frequencies=frequencies

        def _split(arr):
            if arr.ndim == 2 and arr.shape[1] == 2:
                return arr[:, 0], arr[:, 1]
            if arr.ndim == 2 and arr.shape[0] == 2:
                return arr[0, :], arr[1, :]
            raise Exception('spectrum does not have the good shape, should be a 2 dimensional array [velocities, amplitudes]')

        if len(amp) == 0: #if amplitude is not defined its extracted from the spectrum
            self.amp = _split(spectrum)[1]
        else:
            self.amp=np.asarray(amp, dtype=float)
        if len(self.amp) < 2:
            raise ValueError('spectrum {0!r} has {1} channel(s): a Spectrum needs at least 2'.format(
                self.name, len(self.amp)))
        if self.velOrFreq=='vel': #velocity mode
            if len(velocities) == 0:#if the velocities are not defined they are extracted from the spectrum
                if len(spectrum) > 0:
                    self.velocities = _split(spectrum)[0]
                else:   #if the spectrum is not defined either (so if only the amplitude are),
                        #the velocities are defined as a range centered on zero, of size=len(amp) and binsize=1
                    # (half up, as Python 2's round: odd lengths keep all channels)
                    self.velocities=np.arange(-int(len(self.amp)/2),LineStacker.round_half_up(len(self.amp)/2.))
        elif self.velOrFreq=='freq': #frequency mode
            if len(frequencies) == 0: #if the frequencies are not defined they are extracted from the spectrum
                if len(spectrum) > 0:
                    self.frequencies = _split(spectrum)[0]
                else: #if the spectrum is not defined either (so if only the amplitude are),
                        #the frequencies are defined as a range sarting from 0, of size=len(amp) and binsize=1
                    self.frequencies=np.arange(0,len(self.amp))
        axis = self.velocities if self.velOrFreq == 'vel' else self.frequencies
        what = 'velocities' if self.velOrFreq == 'vel' else 'frequencies'
        if len(axis) != len(self.amp):
            raise ValueError('spectrum {0!r}: {1} fluxes but {2} {3}'.format(
                self.name, len(self.amp), len(axis), what))
        try:
            axis = np.asarray(axis, dtype=float)
        except (TypeError, ValueError):
            raise ValueError('spectrum {0!r}: its {1} must be numbers'.format(self.name, what))
        if axis.ndim != 1 or len(axis) < 2:
            raise ValueError('spectrum {0!r}: needs at least 2 {1} (got {2})'.format(
                self.name, what, axis.size))
        if not np.all(np.isfinite(axis)):
            raise ValueError('spectrum {0!r}: its {1} must all be finite'.format(self.name, what))
        steps = np.diff(axis)
        if not (np.all(steps > 0) or np.all(steps < 0)):
            raise ValueError('spectrum {0!r}: its {1} must increase or decrease steadily (no repeated '
                             'values)'.format(self.name, what))
        # channels are located with one channel width: the axis must be evenly
        # spaced (rounding of the values in text files is tolerated)
        step = (axis[-1] - axis[0]) / (len(axis) - 1.)
        offset = np.abs(axis - (axis[0] + step * np.arange(len(axis)))) / abs(step)
        if offset.max() > MAX_AXIS_DEVIATION:
            k = int(np.argmax(offset))
            raise ValueError('spectrum {0!r}: its {1} must be evenly spaced (channel {2} is {3:.2g} channel '
                             'away from an even grid); resample the spectrum to an even grid first'.format(
                                 self.name, what, k, offset[k]))
        if self.velOrFreq == 'vel':
            self.velocities = axis
            self.velBin = step                       # velocity channel width
        else:
            self.frequencies = axis
            self.freqBin = step                      # frequency channel width

        if self.z is not None: #use z and rest emission frequency to find stack center
            if self.fEmLine is False or self.fEmLine is None:
                raise Exception('No emission frequency for spectrum '+repr(self.name)+', you need an emission frequency (fEmLine) when using z')
            fObsLine=self.fEmLine/(1.+self.z)
            if velOrFreq=='vel':
                # velocities (km/s) relative to the rest frequency fEmLine
                vLine = self.freqToVel(fObsLine, 0., self.fEmLine)
                self.centerIndex = LineStacker.round_half_up((vLine - self.velocities[0]) / self.velBin)
                if len(self.frequencies) == 0:
                    self.frequencies = self.velToFreq(np.asarray(self.velocities, dtype=float), 0., self.fEmLine)
                    self.freqBin = (self.frequencies[-1] - self.frequencies[0]) / (len(self.frequencies) - 1.)
            else:
                self.centerIndex=LineStacker.round_half_up((fObsLine-self.frequencies[0])/self.freqBin)
        elif self.centralVelocity is not None:
            self.centerIndex=LineStacker.round_half_up((self.centralVelocity-self.velocities[0])/self.velBin)
        elif self.centralFrequency is not None:
            self.centerIndex=LineStacker.round_half_up((self.centralFrequency-self.frequencies[0])/self.freqBin)
        if fit:#if fit is on, line parameters are estimated using gaussian fitting
            self.fitAgain()

        if isinstance(weights, str) and weights=='1/A': #weights set to 1 over amplitude
            if fit==False: #/!\ if not fitted the amplitude is defined as maximum value /!\
                lineAmp = np.nanmax(self.amp) if np.any(np.isfinite(self.amp)) else np.nan
            else:
                lineAmp = self.gaussfit[1][0]
            if not np.isfinite(lineAmp) or lineAmp <= 0:
                import warnings
                warnings.warn("spectrum {0!r}: weights='1/A' but the line amplitude is {1}; weight set to 0".format(
                    self.name, lineAmp), stacklevel=2)
                self.weights = 0.
            else:
                self.weights = 1. / lineAmp
        elif isinstance(weights, str) and weights=='sigma2': #weights set to 1 over sigma**2
            if fit:
                # lineWidth > 0; the channel width is negative for decreasing axes
                if velOrFreq=='vel':
                    lineWidthInBins=self.lineWidth/abs(self.velBin)
                else:
                    lineWidthInBins=self.lineWidth/abs(self.freqBin)
                leftLimSpectra=int(self.centerIndex-lineWidthInBins*2.35)
                rightLimSpectra=int(self.centerIndex+lineWidthInBins*2.35)
                toSTD=[]
                if rightLimSpectra<len(self.amp)-1:
                    toSTD.append(self.amp[rightLimSpectra:])
                if leftLimSpectra>0:
                    toSTD.append(self.amp[:leftLimSpectra])
                if len(toSTD) == 0:
                    raise Exception('it seems your line is too large to use sigma2 weighting')
                self.weights=1./np.nanstd(np.concatenate(toSTD))**2
            else:
                self.weights=1./np.nanstd(self.amp)**2 #calculating standard deviation of the ENTIRE spectrum
        else:
            self.weights=weights

    #
    def velToFreq(self, vel, z, fEmLine):
        """
            Frequency of the velocity vel (km/s) relative to the frequency
            fEmLine/(1+z) (z=0: the rest frequency), in the spectrum's
            velConvention. Same unit as fEmLine.
        """
        return velocity_to_frequency(vel, fEmLine / (1. + z), self.velConvention)

    def freqToVel(self, freq, z, fEmLine):
        """
            Velocity (km/s) of the frequency freq relative to fEmLine/(1+z)
            (z=0: the rest frequency), in the spectrum's velConvention.
        """
        return frequency_to_velocity(freq, fEmLine / (1. + z), self.velConvention)

    def fitAgain(self):
        if self.velOrFreq=='vel':
            self.gaussfit=LineStacker.tools.fit.GaussFit(     fctToFit=self.amp,
                                                fullFreq=self.velocities,
                                                returnInfos=True)
            self.velCenter=self.gaussfit[1][1]
            self.centerIndex=LineStacker.round_half_up((self.velCenter-self.velocities[0])/self.velBin)
        if self.velOrFreq=='freq':
            self.gaussfit=LineStacker.tools.fit.GaussFit(     fctToFit=self.amp,
                                                fullFreq=self.frequencies,
                                                returnInfos=True)
            self.velCenter=self.gaussfit[1][1]
            self.centerIndex=LineStacker.round_half_up((self.velCenter-self.frequencies[0])/self.freqBin)
        self.lineWidth=abs(self.gaussfit[1][2])


def velocity_to_frequency(vel, fRef, convention='radio'):
    """Frequency (unit of fRef) of velocities vel (km/s) relative to fRef:
    radio f = fRef (1 - v/c), optical f = fRef / (1 + v/c), relativistic
    f = fRef sqrt((1 - v/c) / (1 + v/c))."""
    b = np.asarray(vel, dtype=float) / c
    if convention == 'radio':
        return fRef * (1. - b)
    if convention == 'optical':
        return fRef / (1. + b)
    if convention == 'relativistic':
        return fRef * np.sqrt((1. - b) / (1. + b))
    raise ValueError('convention must be one of {0}'.format(VELOCITY_CONVENTIONS))


def frequency_to_velocity(freq, fRef, convention='radio'):
    """Velocity (km/s) of frequencies freq relative to fRef (inverse of
    velocity_to_frequency)."""
    r = np.asarray(freq, dtype=float) / fRef
    if convention == 'radio':
        return c * (1. - r)
    if convention == 'optical':
        return c * (1. / r - 1.)
    if convention == 'relativistic':
        return c * (1. - r**2) / (1. + r**2)
    raise ValueError('convention must be one of {0}'.format(VELOCITY_CONVENTIONS))


def _set_centers(spectra, center='lineCenterIndex', regridFromZ=False, regridMethod='scaleToMin'):
    """Fix spectrum.centerIndex for every spectrum according to the centring method."""
    for (i,spectrum) in enumerate(spectra):
        if isinstance(center, str) and center=='fit':
            axis = spectrum.velocities if spectrum.velOrFreq == 'vel' else spectrum.frequencies
            if len(axis) > 0:
                step = spectrum.velBin if spectrum.velOrFreq == 'vel' else spectrum.freqBin
                tempGaussfit=LineStacker.tools.fit.GaussFit(     fctToFit=spectrum.amp,
                                                    fullFreq=axis,
                                                    returnInfos=True)[1]
                spectrum.centerIndex=LineStacker.round_half_up((tempGaussfit[1]-axis[0])/step)
            else:
                tempGaussfit=LineStacker.tools.fit.GaussFit(     fctToFit=spectrum.amp,
                                                    returnInfos=True)[1]
                spectrum.centerIndex=LineStacker.round_half_up(tempGaussfit[1])

        elif isinstance(center, str) and center=='zero_vel': #center on v=0km/s
            vel = np.asarray(spectrum.velocities, dtype=float)
            if len(vel) == 0:
                raise ValueError("center='zero_vel' needs velocities (spectrum {0!r} has none)".format(spectrum.name))
            hits = np.where(vel == 0.)[0]
            if len(hits):
                spectrum.centerIndex = int(hits[0])
            else: #if the value 0 is not found, search closest
                centerIndex=-vel[0]/spectrum.velBin
                if centerIndex < -0.5 or centerIndex > len(vel) - 0.5:
                    raise ValueError('zero velocity is not in the velocities of spectrum {0!r}'.format(spectrum.name))
                spectrum.centerIndex=LineStacker.round_half_up(centerIndex)

        elif isinstance(center, str) and center=='center': #center on the center of the array
            spectrum.centerIndex=int(len(spectrum.amp)/2)

        elif isinstance(center, str) and center=='lineCenterIndex':
            if spectrum.centerIndex is None:
                raise Exception('center set to lineCenterIndex but no center index is found for line '+str(i))
        else: #user defined center
            spectrum.centerIndex=int(center)
            if regridFromZ:
                if regridMethod=='scaleToMin':
                    zToScale=min([tempSp.z for tempSp in spectra])
                elif regridMethod=='scaleToMax':
                    zToScale=max([tempSp.z for tempSp in spectra])
                zRatio=(1.+spectrum.z)/(1.+zToScale)
                spectrum.centerIndex=LineStacker.round_half_up(center*zRatio)
        spectrum.centerIndex = int(spectrum.centerIndex)


class Aligned(object):
    """All spectra aligned on their line centre.

    X       (nspectra, L) fluxes, NaN where a spectrum has no channel
    W       (nspectra, L) weights, 0 where X is NaN
    lo, hi  per spectrum: channels [lo, hi) of the full alignment are its own
            'full' stack range (as computed by Stack for that spectrum alone)
    """

    def __init__(self, spectra, chansStack='full'):
        n = len(spectra)
        ci = np.array([sp.centerIndex for sp in spectra], dtype=int)
        namp = np.array([len(sp.amp) for sp in spectra], dtype=int)
        out = np.where((ci < 0) | (ci >= namp))[0]
        if out.size:
            i = int(out[0])
            raise ValueError('spectrum {0} ({1!r}): its line centre, channel {2}, is outside the spectrum '
                             '(channels 0 to {3}); check its z, central velocity/frequency or centerIndex'.format(
                                 i, spectra[i].name, ci[i], namp[i] - 1))
        # Stack sizes its 'full' output with the length of the spectral axis
        # of the last spectrum's kind (velocities, else frequencies).
        last = spectra[-1]
        if len(last.velocities) > 0:
            nspec = np.array([len(sp.velocities) for sp in spectra], dtype=int)
        else:
            nspec = np.array([len(sp.frequencies) for sp in spectra], dtype=int)
        if chansStack == 'full':
            left, right = int(ci.max()), int((nspec - ci).max())
            self.offset0 = -left
            L = left + right
        else:
            L = int(chansStack)
            self.offset0 = -int(L / 2.)
        self.L = L
        self.full = chansStack == 'full'
        self.X = np.full((n, L), np.nan)
        self.W = np.zeros((n, L))
        for i, sp in enumerate(spectra):
            src = ci[i] + self.offset0 + np.arange(L)          # channel of spectrum i
            ok = (src >= 0) & (src < namp[i])
            amp = np.asarray(sp.amp, dtype=float)
            self.X[i, ok] = amp[src[ok]]
            w = sp.weights
            if np.ndim(w) == 0:
                self.W[i, ok] = float(w)
            else:
                w = np.asarray(w, dtype=float)
                okw = ok & (src < len(w))
                self.W[i, okw] = w[src[okw]]
        self.M = np.isfinite(self.X)
        self.X[~self.M] = np.nan               # inf counts as missing, for medians too
        # channels without a finite flux (NaN, inf) carry no weight
        self.W[~self.M] = 0.
        bad = ~np.isfinite(self.W) | (self.W < 0)
        if bad.any():
            i = int(np.where(bad.any(axis=1))[0][0])
            raise ValueError('spectrum {0} ({1!r}) has weight {2}: weights must be finite and >= 0'.format(
                i, spectra[i].name, self.W[i][bad[i]][0]))
        self.X0 = np.where(self.M, self.X, 0.)
        # each spectrum's own 'full' range within the alignment
        self.lo = (-ci) - self.offset0
        self.hi = (nspec - ci) - self.offset0


def _median_sorted(values, valid_count):
    """Median along axis -2 of values whose NaNs are sorted last (NaN where
    there is no valid value)."""
    k = valid_count
    lo = np.clip((k - 1) // 2, 0, None)
    hi = np.clip(k // 2, 0, None)
    a = np.take_along_axis(values, lo[..., None, :], axis=-2)[..., 0, :]
    b = np.take_along_axis(values, hi[..., None, :], axis=-2)[..., 0, :]
    out = 0.5 * (a + b)
    return np.where(k > 0, out, np.nan)


def combine_counts(al, counts, method='mean'):
    """Stacks of many resamplings at once.

    counts : (nsets, nspectra) how many times each spectrum enters each set
    returns (stacks (nsets, L), sources-per-channel (nsets, L))
    """
    counts = np.asarray(counts, dtype=float)
    nsrc = counts @ al.M
    if method == 'mean':
        num = counts @ (al.W * al.X0)
        den = counts @ al.W
        with np.errstate(invalid='ignore', divide='ignore'):
            stacks = np.where(den > 0, num / den, np.nan)       # NaN: no data
        return stacks, nsrc
    if method == 'median':
        raise ValueError('use combine_indices for medians')
    raise ValueError("method must be 'mean' or 'median'")


def combine_indices(al, idx, method='mean'):
    """Stacks of resamplings given as index arrays (nsets, k)."""
    idx = np.asarray(idx)
    if method == 'mean':
        counts = np.zeros((idx.shape[0], al.X.shape[0]))
        np.add.at(counts, (np.repeat(np.arange(idx.shape[0]), idx.shape[1]), idx.ravel()), 1.)
        return combine_counts(al, counts, 'mean')
    if method == 'median':
        sel = al.X[idx]                                   # (nsets, k, L)
        k = (~np.isnan(sel)).sum(axis=1)
        sel = np.sort(sel, axis=1)                        # NaNs last
        return _median_sorted(sel, k), k.astype(float)
    raise ValueError("method must be 'mean' or 'median'")


def crop(al, stack, members):
    """Cut a full-alignment stack to the range Stack() returns for the
    spectra 'members' (no-op for a fixed chansStack)."""
    if not al.full:
        return stack
    lo = int(al.lo[members].min())
    hi = int(al.hi[members].max())
    return stack[..., lo:hi]


def Stack(  spectra,
            chansStack='full',
            method='mean',
            center='lineCenterIndex',
            regridFromZ=False,
            regridMethod='scaleToMin',
            verbose=None):

    """
        Stacks spectra (LineStacker.OneD_Stacker.Spectrum objects), aligned on
        their lines: weighted mean or median, channel by channel.

        Parameters
        ---------
        spectra
            List of spectra (LineStacker.OneD_Stacker.Spectrum objects).
        chansStack
            Number of channels to stack. Set to 'full' to stack all channels from all spectra. User input (int) otherwise
        method
            stacking method, **'mean'** and **'median'** supported
        center
            Method to find central frequency of the stack, possible values are:\n
            **"center"**, to stack all spectra center to center,\n
            **'fit'** to use gaussian fitting on the spectrum to determine line center,\n
            **'zero_vel'** to stack on velocity=0 bin,\n
            **'lineCenterIndex'** use the line center initiated with the spectrum,\n
            Or directly defined by the user (int)
        regridFromZ
            if set to True spectra will be regridded to take into account the redshift difference of the different sources. See LineStacker.analysisTools.regridFromZ1D for a more complete description.
        regridMethod
            Used if regridFromZ is True. Can be set either to **'scaleToMin'** or **'scaleToMax'**.
            In the first case all spectra are regridded to match the smallest redshift (over-gridding), all spectra are regridded to match the highest redshift in the other case (under-gridding).
        verbose
            report the regridding (default: LineStacker.verbose)

        Returns
        -------
        tuple
            (stacked spectrum, number of spectra in each channel). Channels
            without valid data (no spectrum covers them, or only blank
            fluxes or zero weights) are NaN.
    """

    if chansStack!='full' and not isinstance(chansStack, (int, np.integer)):
        raise Exception("chansStack must be an int or 'full' for full spectra" )
    spectra = list(spectra)
    if not spectra:
        raise ValueError('no spectra to stack')
    if method not in ('mean', 'median'):
        raise ValueError("method must be 'mean' or 'median'")
    if regridFromZ:
        import LineStacker.analysisTools
        spectra=LineStacker.analysisTools.regridFromZ1D(spectra, regridMethod=regridMethod, verbose=verbose)
    _set_centers(spectra, center, regridFromZ, regridMethod)
    al = Aligned(spectra, chansStack)
    everyone = np.arange(len(spectra))[None, :]
    stacks, nsrc = combine_indices(al, everyone, method)
    return stacks[0], nsrc[0]
