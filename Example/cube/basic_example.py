import LineStacker
import LineStacker.line_image

coords=LineStacker.readCoordsGUI()
imagenames=LineStacker.selectImagesGUI()

stacked=LineStacker.line_image.stack(   coords,
                                        imagenames=imagenames)
