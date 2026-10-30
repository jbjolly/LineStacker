import sys
import glob

import LineStacker
import LineStacker.image


def main():
    #coordinates files are selected using the GUI,
    #or all *_coord.txt files of a folder given on the command line
    if len(sys.argv) > 1:
        coordNames=sorted(glob.glob(sys.argv[1].rstrip('/')+'/*_coord.txt'))
    else:
        coordNames=LineStacker.readCoordsNamesGUI()
    coords=LineStacker.readCoords(coordNames)

    #image names are identical to coordinates files,
    #with '.image' replacing '_coord.txt'
    imagenames=([coord.removesuffix('_coord.txt')+'.image' for coord in coordNames])

    #because redshift is used to identify the line center,
    #the emission frequency is also provided
    stacked=LineStacker.image.stack(   coords,
                                            imagenames=imagenames,
                                            fEm=1897420620253.1646,
                                            stampsize=16,
                                            N_chans=16)
    return stacked


if __name__ == '__main__':
    main()
