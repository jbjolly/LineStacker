import LineStacker
import LineStacker.line_image

coordNames=LineStacker.readCoordsNamesGUI()
coords=LineStacker.readCoords(coordNames)
imagenames=([coord.strip('_coords.txt')+'.image' for coord in coordNames])

stacked=LineStacker.line_image.stack(   coords,
                                        imagenames=imagenames,
                                        fEm=1900.548e9,
                                        stampsize=16)
