import LineStacker
import LineStacker.line_image

coordNames=LineStacker.readCoordsNamesGUI()
coords=LineStacker.readCoords(coordNames)
imagenames=([coord.strip('_coords.txt')+'.image' for coord in coordNames])

stacked=LineStacker.line_image.stack(   coords,
                                        imagenames=imagenames,
                                        fEm=1900.548e9,
                                        stampsize=16)

import matplotlib.pyplot as plt

fig=plt.figure()
ax=fig.add_subplot(111)
ax.plot(stacked[0])
fig.show()
