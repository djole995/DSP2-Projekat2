
#include "ImageProcessing.h"
#include "ImageInterpolation.h"

#include <QDebug>

void imageProcessingFun(const QString& progName, QImage* const outImgs, const QImage* const inImgs, const QVector<double>& params) 
{
	int X_SIZE = inImgs->width();
	int Y_SIZE = inImgs->height();

	/* NOTE: Calculate output image resolution and construct output image object */

	if(progName == "Sample and hold") 
	{	
		/* Input image data in RGB format can be obtained with inImgs->bits() */
		/* Vertical scale factor is params[0] */
		/* Horizontal scale factor is params[1] */

		int newXSize = ((int) (X_SIZE * params[1])) % 4 == 0 ? (int)(X_SIZE * params[1]) : (int)(X_SIZE * params[1]) + 4-((int)(X_SIZE * params[1])) % 4;
		int newYSize = ((int)(Y_SIZE * params[0])) % 4 == 0 ? (int)(Y_SIZE * params[0]) : (int)(Y_SIZE * params[0]) + 4 - ((int)(Y_SIZE * params[0])) % 4;

		/* TO DO: Calculate output image resolution and construct output image object */
		*outImgs = *(new QImage(newXSize, newYSize, inImgs->format()));
		/* TO DO: Perform Sample and hold interpolation  */
		sampleAndHold(inImgs->bits(), X_SIZE, Y_SIZE, outImgs->bits(), newXSize,  newYSize);
	}
	else if (progName == "Bilinear") 
	{
		/* Input image data in RGB format can be obtained with inImgs->bits() */
		/* Vertical scale factor is params[0] */
		/* Horizontal scale factor is params[1] */
		int newXSize = ((int)(X_SIZE * params[1])) % 4 == 0 ? (int)(X_SIZE * params[1]) : (int)(X_SIZE * params[1]) + 4 - ((int)(X_SIZE * params[1])) % 4;
		int newYSize = ((int)(Y_SIZE * params[0])) % 4 == 0 ? (int)(Y_SIZE * params[0]) : (int)(Y_SIZE * params[0]) + 4 - ((int)(Y_SIZE * params[0])) % 4;
		/* TO DO: Calculate output image resolution and construct output image object */
		*outImgs = *(new QImage(newXSize, newYSize, inImgs->format()));
		bilinearInterpolate(inImgs->bits(), X_SIZE, Y_SIZE, outImgs->bits(), newXSize, newYSize);
	}
	else if (progName == "Bicubic")
	{
		/* Input image data in RGB format can be obtained with inImgs->bits() */
		/* Vertical scale factor is params[0] */
		/* Horizontal scale factor is params[1] */

		/* TO DO: Calculate output image resolution and construct output image object */

		/* TO DO: Perform Bicubic interpolation  */

	}
	else if(progName == "Rotation") 
	{	
		/* Input image data in RGB format can be obtained with inImgs->bits() */
		/* Rotation angle in degrees is params[0]*/
		/* Center of rotation coordinates are (XSIZE/2, YSIZE/2) */

		/* TO DO: Construct output image object */

		/* TO DO: Perform image rotation */
	
	}
	else if (progName == "Rotation Bilinear") 
	{
		/* Input image data in RGB format can be obtained with inImgs->bits() */
		/* Rotation angle in degrees is params[0]*/
		/* Center of rotation coordinates are (XSIZE/2, YSIZE/2) */

		/* TO DO: Construct output image object */

		/* TO DO: Perform image rotation with bilinear interpolation */
	}

}

