#include "ImageInterpolation.h"
#include "ColorSpaces.h"
#include <math.h>


void sampleAndHold(const uchar input[], int xSize, int ySize, uchar output[], int newXSize, int newYSize)
{
	double xFactor = (double) newXSize / xSize;
	double yFactor = (double) newYSize / ySize;
	int ii, jj;
	double yPos;
	double xPos;

	uchar *Y_buff_in = new uchar[xSize*ySize];
	char *U_buff_in = new char[xSize/2*ySize/2];
	char *V_buff_in = new char[xSize / 2 * ySize / 2];

	uchar *Y_buff_out = new uchar[newXSize*newYSize];
	char *U_buff_out = new char[newXSize / 2 * newYSize / 2];
	char *V_buff_out = new char[newXSize / 2 * newYSize / 2];

	RGBtoYUV420(input, xSize, ySize, Y_buff_in, U_buff_in, V_buff_in);

	for (int i = 0; i < newYSize; i++)
	{
		yPos = i / yFactor;

		if (yPos - floor(yPos) < 0.5)
			ii = (int)floor(yPos);
		else
			ii = (int)ceil(yPos);

		for (int j = 0; j < newXSize; j++)
		{
			xPos = j / xFactor;

			if (xPos - floor(xPos) < 0.5)
				jj = (int)floor(xPos);
			else
				jj = (int)ceil(xPos);

			Y_buff_out[i*newXSize + j] = Y_buff_in[ii*xSize + jj];
		}
	}

	for (int i = 0; i < newYSize/2; i++)
	{
		yPos = i / yFactor;

		if (yPos - floor(yPos) < 0.5)
			ii = (int)floor(yPos);
		else
			ii = (int)ceil(yPos);

		for (int j = 0; j < newXSize/2; j++)
		{
			xPos = j / xFactor;

			if (xPos - floor(xPos) < 0.5)
				jj = (int)floor(xPos);
			else
				jj = (int)ceil(xPos);

			U_buff_out[i*newXSize/2 + j] = U_buff_in[ii*xSize/2 + jj];
			V_buff_out[i*newXSize/2 + j] = V_buff_in[ii*xSize/2 + jj];
		}
	}

	YUV420toRGB(Y_buff_out, U_buff_out, V_buff_out, newXSize, newYSize, output);
}

void bilinearInterpolate(const uchar input[], int xSize, int ySize, uchar output[], int newXSize, int newYSize)
{
	double xFactor = (double)newXSize / xSize;
	double yFactor = (double)newYSize / ySize;
	int ii, jj;
	double yPos;
	double xPos;

	double a, b;

	uchar *Y_buff_in = new uchar[xSize*ySize];
	char *U_buff_in = new char[xSize / 2 * ySize / 2];
	char *V_buff_in = new char[xSize / 2 * ySize / 2];

	uchar *Y_buff_out = new uchar[newXSize*newYSize];
	char *U_buff_out = new char[newXSize / 2 * newYSize / 2];
	char *V_buff_out = new char[newXSize / 2 * newYSize / 2];

	RGBtoYUV420(input, xSize, ySize, Y_buff_in, U_buff_in, V_buff_in);

	for (int i = 0; i < newYSize; i++)
	{
		yPos = i / yFactor;

		ii = (int) floor(yPos);

		b = i / yFactor - floor(i / yFactor);

		for (int j = 0; j < newXSize; j++)
		{
			a = j / xFactor - floor(j / xFactor);
			xPos = j / xFactor; 

			jj = (int)floor(xPos);

			Y_buff_out[i*newXSize + j] = Y_buff_in[ii*xSize + jj]*(1-a)*(1-b) + 
										 Y_buff_in[(ii+1)*xSize + jj] * (1 - a)*b
										 + Y_buff_in[ii*xSize + jj+1] * a*(1 - b)
										 + Y_buff_in[(ii+1)*xSize + jj + 1] * a*b;
		}
	}

	for (int i = 0; i < newYSize/2; i++)
	{
		yPos = i / yFactor;

		ii = (int)floor(yPos);

		for (int j = 0; j < newXSize/2; j++)
		{
			xPos = j / xFactor;

			jj = (int)floor(xPos);

			U_buff_out[i*newXSize / 2 + j] = U_buff_in[ii*xSize / 2 + jj] * (1 - a)*(1 - b)+
				U_buff_in[(ii + 1)*xSize/2 + jj] *(1 - a)*b
				+ U_buff_in[ii*xSize / 2 + jj + 1] *a*(1 - b);
				+ U_buff_in[(ii + 1)*xSize/2 + jj + 1] *a*b;

			V_buff_out[i*newXSize / 2 + j] = V_buff_in[ii*xSize / 2 + jj] * (1 - a)*(1 - b)+
				V_buff_in[(ii + 1)*xSize/2 + jj] * (1 - a)*b
				+ V_buff_in[ii*xSize / 2 + jj + 1] * a*(1 - b);
				+ V_buff_in[(ii + 1)*xSize/2 + jj + 1] *a*b;
		}
	}

	YUV420toRGB(Y_buff_out, U_buff_out, V_buff_out, newXSize, newYSize, output);
}

void bicubicInterpolate(uchar input[], int xSize, int ySize, uchar output[], int newXSize, int newYSize)
{
	/* TO DO */
}

void imageRotate(const uchar input[], int xSize, int ySize, uchar output[], int m, int n, double angle)
{
	/* TO DO */
}

void imageRotateBilinear(const uchar input[], int xSize, int ySize, uchar output[], int m, int n, double angle)
{
	/* TO DO */
}