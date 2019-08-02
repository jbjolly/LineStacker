import LineStacker
import LineStacker.line_image

#coordinates files are selected using the GUI
coordNames=LineStacker.readCoordsNamesGUI()
coords=LineStacker.readCoords(coordNames)

#image names are identical to coordinates files, with '.image' replacing '_coords.txt'
imagenames=([coord.strip('_coords.txt')+'.image' for coord in coordNames])

#because redshift is used to identify the line center, the emission frequency is also provided
stacked=LineStacker.line_image.stack(   coords,
                                        imagenames=imagenames,
                                        fEm=1897420620253.1646,
                                        stampsize=16)
