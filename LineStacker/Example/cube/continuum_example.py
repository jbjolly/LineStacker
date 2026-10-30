#Continuum stacking with the same function as line stacking.
#The example cubes have a spectral axis, so with continuum=True LineStacker
#warns and stacks the images summed over channels (computed in memory).
#Images without a spectral axis are stacked as continuum automatically.
import sys
import glob

import LineStacker
import LineStacker.image
import LineStacker.analysisTools


def main():
    if len(sys.argv) > 1:
        coordNames=sorted(glob.glob(sys.argv[1].rstrip('/')+'/*_coord.txt'))
    else:
        coordNames=LineStacker.readCoordsNamesGUI()
    coords=LineStacker.readCoords(coordNames)
    imagenames=([coord.removesuffix('_coord.txt')+'.image' for coord in coordNames])

    stacked=LineStacker.image.stack(    coords,
                                        imagenames=imagenames,
                                        stampsize=16,
                                        continuum=True,
                                        outfile='continuumStack.fits')
    print('stacked continuum flux at the centre: {0:.3g}'.format(stacked[0]))

    #how often would random positions give as much? (all CPU cores)
    random=LineStacker.analysisTools.stack_estimator(  coords,
                                                        imagenames=imagenames,
                                                        continuum=True,
                                                        nRandom=1000,
                                                        lowerLimit=4e-6,
                                                        upperLimit=9e-6,
                                                        n_jobs=-1)
    print('random positions: mean {0:.3g}, std {1:.3g}'.format(random[1], random[0]))
    return stacked, random


if __name__ == '__main__':
    main()
