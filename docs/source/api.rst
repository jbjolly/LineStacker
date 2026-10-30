API reference
=============

LineStacker
-----------

.. automodule:: LineStacker
   :members: Coord, readCatalog, makeCoords, readCoords, readCoordsNamesGUI, getPixelCoords, randomCoords, randomizeCoords, round_half_up, verbose, casa_logfile, cache_check

LineStacker.image
-----------------

.. automodule:: LineStacker.image
   :members: stack, GaussianPrimaryBeam, readCube, getFlux, clear_cache

LineStacker.OneD_Stacker
------------------------

.. automodule:: LineStacker.OneD_Stacker
   :members: Spectrum, Stack, velocity_to_frequency, frequency_to_velocity

LineStacker.analysisTools
-------------------------

.. automodule:: LineStacker.analysisTools
   :members: bootstraping_OneD, bootStraping_Cube, subsample_OneD, maximizeSNR, maximizeAmp, maximizeOutflow, stack_estimator, noise_estimator, randomizeCoords, rebin_OneD, rebin_CubesSpectra, regridFromZ, regridFromZ1D, congrid

LineStacker.tools
-----------------

.. automodule:: LineStacker.tools
   :members: ProgressBar, resolve_n_jobs

LineStacker.tools.fit
---------------------

.. automodule:: LineStacker.tools.fit
   :members: GaussFit, DoubleGaussFit, DoublePeakFit, FitError, gaussFct, Gauss, DoubleGauss, DoubleGaussPeak

LineStacker.tools.parallel
--------------------------

.. automodule:: LineStacker.tools.parallel
   :members: map_chunks, get_shared, set_default_backend
