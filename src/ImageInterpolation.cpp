#include "ImageInterpolation.h"
#include "ColorSpaces.h"
#include "ImageFilter.h"
#include <math.h>

void sampleAndHold(const uchar input[], int xSize, int ySize, uchar output[], int newXSize, int newYSize)
{
	double xFactor = (double) newXSize / xSize;
	double yFactor = (double) newYSize / ySize;
	int ii, jj;

	uchar *Y_buff_in = new uchar[xSize*ySize];
	char *U_buff_in = new char[xSize/2*ySize/2];
	char *V_buff_in = new char[xSize / 2 * ySize / 2];

	uchar *Y_buff_out = new uchar[newXSize*newYSize];
	char *U_buff_out = new char[newXSize / 2 * newYSize / 2];
	char *V_buff_out = new char[newXSize / 2 * newYSize / 2];

	RGBtoYUV420(input, xSize, ySize, Y_buff_in, U_buff_in, V_buff_in);

	for (int i = 0; i < newYSize; i++)
	{
		ii = (int) floor((i-1) / yFactor + 1 + 0.5);

		if (ii < 0)
			ii = 0;
		else if (ii > ySize - 1)
			ii = ySize - 1;

		for (int j = 0; j < newXSize; j++)
		{
			jj = (int) floor((j - 1) / xFactor + 1 + 0.5);

			if (jj < 0)
				jj = 0;
			else if (jj > xSize - 1)
				jj = xSize - 1;

			Y_buff_out[i*newXSize + j] = Y_buff_in[ii*xSize + jj];
		}
	}

	for (int i = 0; i < newYSize/2; i++)
	{
		ii = (int)floor((i - 1) / yFactor + 1 + 0.5);

		if (ii < 0)
			ii = 0;
		else if (ii > ySize/2 - 1)
			ii = ySize/2 - 1;

		for (int j = 0; j < newXSize/2; j++)
		{
			jj = (int)floor((j - 1) / xFactor + 1 + 0.5);

			if (jj < 0)
				jj = 0;
			else if (jj > xSize/2 - 1)
				jj = xSize/2 - 1;

			U_buff_out[i*newXSize/2 + j] = U_buff_in[ii*xSize/2 + jj];
			V_buff_out[i*newXSize/2 + j] = V_buff_in[ii*xSize/2 + jj];
		}
	}

	YUV420toRGB(Y_buff_out, U_buff_out, V_buff_out, newXSize, newYSize, output);

	delete Y_buff_in;
	delete U_buff_in;
	delete V_buff_in;
	delete Y_buff_out;
	delete U_buff_out;
	delete V_buff_out;
}

void bilinearInterpolate(const uchar input[], int xSize, int ySize, uchar output[], int newXSize, int newYSize)
{
	double xFactor = (double)newXSize / xSize;
	double yFactor = (double)newYSize / ySize;
	int m1, m2, n1, n2;
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
		yPos = (i-1) / yFactor + 1;

		if (yPos > ySize - 1)
			yPos = ySize - 1;

		m1 = (int)floor(yPos);
		m2 = (int)ceil(yPos);

		b = yPos - floor(yPos);

		for (int j = 0; j < newXSize; j++)
		{
			xPos = (j-1) / xFactor + 1;

			if (xPos > xSize - 1)
				xPos = xSize - 1;

			n1 = (int)floor(xPos);
			n2 = (int)ceil(xPos);

			a = xPos - floor(xPos);

			Y_buff_out[i*newXSize + j] = Y_buff_in[m1*xSize + n1]*(1-a)*(1-b) + 
										 Y_buff_in[m2*xSize + n1] * (1 - a)*b
										 + Y_buff_in[m1*xSize + n2] * a*(1 - b)
										 + Y_buff_in[m2*xSize + n2] * a*b;
		}
	}

	for (int i = 0; i < newYSize/2; i++)
	{
		yPos = (i - 1) / yFactor + 1;

		if (yPos > ySize/2 - 1)
			yPos = ySize/2 - 1;

		m1 = (int)floor(yPos);
		m2 = (int)ceil(yPos);

		b = yPos - floor(yPos);

		for (int j = 0; j < newXSize/2; j++)
		{
			xPos = (j-1) / xFactor + 1;

			if (xPos > xSize / 2 - 1)
				xPos = xSize / 2 - 1;

			n1 = (int)floor(xPos);
			n2 = (int)ceil(xPos);

			a = xPos - floor(xPos);

			U_buff_out[i*newXSize / 2 + j] = U_buff_in[m1*xSize / 2 + n1] * (1 - a)*(1 - b)+
				U_buff_in[m2*xSize/2 + n1] *(1 - a)*b
				+ U_buff_in[m1*xSize / 2 + n2] * a*(1 - b)
				+ U_buff_in[m2*xSize/2 + n2] *a*b;

			V_buff_out[i*newXSize / 2 + j] = V_buff_in[m1*xSize / 2 + n1] * (1 - a)*(1 - b) +
				V_buff_in[m2*xSize / 2 + n1] * (1 - a)*b
				+ V_buff_in[m1*xSize / 2 + n2] * a*(1 - b)
				+ V_buff_in[m2*xSize / 2 + n2] * a*b;
		}
	}

	YUV420toRGB(Y_buff_out, U_buff_out, V_buff_out, newXSize, newYSize, output);

	delete Y_buff_in;
	delete U_buff_in;
	delete V_buff_in;
	delete Y_buff_out;
	delete U_buff_out;
	delete V_buff_out;
}

void bicubicInterpolate(const uchar input[], int xSize, int ySize, uchar output[], int newXSize, int newYSize)
{
	double xFactor = (double)newXSize / xSize;
	double yFactor = (double)newYSize / ySize;
	int m[4];
	int n[4];
	double yPos, xPos;
	uchar pointsY[4];
	uchar pointsU[4];
	uchar pointsV[4];
	uchar tmpBuffY[4];
	uchar tmpBuffU[4];
	uchar tmpBuffV[4];

	uchar *Y_buff_in = new uchar[xSize*ySize];
	char *U_buff_in = new char[xSize / 2 * ySize / 2];
	char *V_buff_in = new char[xSize / 2 * ySize / 2];

	uchar *Y_buff_out = new uchar[newXSize*newYSize];
	char *U_buff_out = new char[newXSize / 2 * newYSize / 2];
	char *V_buff_out = new char[newXSize / 2 * newYSize / 2];

	RGBtoYUV420(input, xSize, ySize, Y_buff_in, U_buff_in, V_buff_in);



	for (int i = 0; i < newYSize; i++)
	{
		yPos = (i - 1) / yFactor + 1;

		m[0] = (yPos > 1) ? (int)floor(yPos) - 1 : 0;
		m[1] = (int)floor(yPos);
		m[2] = ((int)yPos == 0) ? 1 : (yPos < ySize-1) ? m[1]+1 : ySize-1;
		m[3] = (yPos < ySize-2) ? m[2] + 1 : ySize-1;

		for (int j = 0; j < newXSize; j++)
		{
			xPos = (j - 1) / xFactor + 1;

			n[0] = (xPos > 1) ? (int)floor(xPos) - 1 : 0;
			n[1] = (int)floor(xPos);
			n[2] = ((int)xPos == 0) ? 1 : (xPos < xSize - 1) ? n[1] + 1 : xSize - 1;
			n[3] = (xPos < xSize - 2) ? n[2] + 1 : xSize - 1;

			for (int ii = 0; ii < 4; ii++)
			{
				for (int jj = 0; jj < 4; jj++)
				{
					pointsY[jj] = Y_buff_in[m[ii] * xSize + n[jj]];
				}
				double d = (xPos > 1) ? xPos - (double)n[0] : xPos - (double)n[0]+1;
				d = (xPos < xSize - 2) ? xPos - (double)n[0] : xPos - (double)(xSize - 3);
				tmpBuffY[ii] = (uchar) cubicInterpolate(pointsY, d, 0);
			}
			
			double d = (yPos > 1) ? yPos - (double)m[0] : yPos - (double)m[0] + 1;
			d = (yPos < ySize - 2) ? yPos - (double)m[0] : yPos - (double)(ySize-3);
			Y_buff_out[i*newXSize + j] = (uchar) cubicInterpolate(tmpBuffY, d, 0);
		}
	}
	
	for (int i = 0; i < newYSize/2; i++)
	{
		yPos = (i - 1) / yFactor + 1;

		m[0] = (yPos > 1) ? (int)floor(yPos) - 1 : 0;
		m[1] = (int)floor(yPos);
		m[2] = ((int)yPos == 0) ? 1 : (yPos < ySize/2 - 1) ? m[1] + 1 : ySize/2 - 1;
		m[3] = (yPos < ySize/2 - 2) ? m[2] + 1 : ySize/2 - 1;

		for (int j = 0; j < newXSize/2; j++)
		{
			xPos = (j - 1) / xFactor + 1;

			n[0] = (xPos > 1) ? (int)floor(xPos) - 1 : 0;
			n[1] = (int)floor(xPos);
			n[2] = ((int)xPos == 0) ? 1 : (xPos < xSize/2 - 1) ? n[1] + 1 : xSize/2 - 1;
			n[3] = (xPos < xSize/2 - 2) ? n[2] + 1 : xSize/2 - 1;

			for (int ii = 0; ii < 4; ii++)
			{
				for (int jj = 0; jj < 4; jj++)
				{
					pointsU[jj] = U_buff_in[m[ii] * xSize / 2 + n[jj]];
					pointsV[jj] = V_buff_in[m[ii] * xSize / 2 + n[jj]];
				}
				double d = (xPos > 1) ? xPos - (double)n[0] : xPos - (double)n[0] + 1;
				d = (xPos < xSize/2 - 2) ? xPos - (double)n[0] : xPos - (double)(xSize/2 - 3);
				tmpBuffU[ii] = (char) cubicInterpolate(pointsU, d, 1);
				tmpBuffV[ii] = (char) cubicInterpolate(pointsV, d, 1);
			}
			double d = (yPos > 1) ? yPos - (double)m[0] : yPos - (double)m[0] + 1;
			d = (yPos < ySize/2 - 2) ? yPos - (double)m[0] : yPos - (double)(ySize/2 - 3);
			U_buff_out[i*newXSize / 2 + j] = (char) cubicInterpolate(tmpBuffU, d, 1);
			V_buff_out[i*newXSize/2 + j] = (char) cubicInterpolate(tmpBuffV, d, 1);
		}
	}

	YUV420toRGB(Y_buff_out, U_buff_out, V_buff_out, newXSize, newYSize, output);

	delete Y_buff_in;
	delete U_buff_in;
	delete V_buff_in;
	delete Y_buff_out;
	delete U_buff_out;
	delete V_buff_out;
}

int cubicInterpolate(uchar points[], double d, int type)
{
	double w[4];
	double distance[3];
	double a;
	int ret = 0;

	a = abs(d);

	if (a < 1)
		w[0] = 3.0 / 2 * a*a*a - 5.0 / 2 * a*a + 1;
	else if (a < 2)
		w[0] = -1.0 / 2 * a*a*a + 5.0 / 2 * a*a - 4 * a + 2;
	else
		w[0] = 0;

	for (int i = 1; i < 4; i++)
	{
		distance[i-1] = abs(i - d);

		a = distance[i - 1];

		if (a < 1)
			w[i] = 3.0 / 2 * a*a*a - 5.0 / 2 * a*a + 1;
		else if (a < 2)
			w[i] = -1.0 / 2 * a*a*a + 5.0 / 2 * a*a - 4 * a + 2;
		else
			w[i] = 0;
	}

	/* Y processing (unsigned char). */
	if (type == 0)
	{
		for (int i = 0; i < 4; i++)
		{
			ret += w[i] * points[i];
		}

		if (ret > 255)
			ret = 255;
		else if (ret < 0)
			ret = 0;
	}

	/* UV processing (signed char). */
	else
	{
		for (int i = 0; i < 4; i++)
		{
			ret += w[i] * (char)points[i];
		}

		if (ret > 127)
			ret = 127;
		else if (ret < -128)
			ret = -128;
	}





	return ret;
}

void imageRotate(const uchar input[], int xSize, int ySize, uchar output[], int m, int n, double angle)
{
	int ii, jj;

	uchar *Y_buff_in = new uchar[xSize*ySize];
	char *U_buff_in = new char[xSize / 2 * ySize / 2];
	char *V_buff_in = new char[xSize / 2 * ySize / 2];

	uchar *Y_buff_out = new uchar[xSize*ySize];
	char *U_buff_out = new char[xSize / 2 * ySize / 2];
	char *V_buff_out = new char[xSize / 2 * ySize / 2];

	RGBtoYUV420(input, xSize, ySize, Y_buff_in, U_buff_in, V_buff_in);



	for (int i = 0; i < ySize; i++)
	{
		for (int j = 0; j < xSize; j++)
		{
			jj = floor(j*cos(angle) - i*sin(angle) + m*sin(angle) - n*cos(angle) + n + 0.5);
			ii = floor(i*cos(angle) + j*sin(angle) - m*cos(angle) - n*sin(angle) + m + 0.5);

			if (ii > ySize - 1 || ii < 0 || jj > xSize-1 || jj < 0)
				Y_buff_out[i*xSize + j] = 0;
			else
				Y_buff_out[i*xSize + j] = Y_buff_in[ii*xSize+jj];
		}
	}

	for (int i = 0; i < ySize/2; i++)
	{
		for (int j = 0; j < xSize/2; j++)
		{
			jj = floor(j*cos(angle) - i*sin(angle) + m/2.0*sin(angle) - n/2.0*cos(angle) + n/2.0 + 0.5);
			ii = floor(i*cos(angle) + j*sin(angle) - m/2.0*cos(angle) - n/2.0*sin(angle) + m/2.0 + 0.5);

			if (ii > ySize/2 - 1 || ii < 0 || jj > xSize/2 - 1 || jj < 0)
			{
				U_buff_out[i*xSize/2 + j] = 0;
				V_buff_out[i*xSize/2 + j] = 0;
			}
			else
			{
				U_buff_out[i*xSize/2 + j] = U_buff_in[ii*xSize/2 + jj];
				V_buff_out[i*xSize/2 + j] = V_buff_in[ii*xSize/2 + jj];
			}
		}
	}

	YUV420toRGB(Y_buff_out, U_buff_out, V_buff_out, xSize, ySize, output);

	delete Y_buff_in;
	delete U_buff_in;
	delete V_buff_in;
	delete Y_buff_out;
	delete U_buff_out;
	delete V_buff_out;
}

void imageRotateBilinear(const uchar input[], int xSize, int ySize, uchar output[], int m, int n, double angle)
{
	float ii, jj;
	int m1, m2, n1, n2;
	double yPos;
	double xPos;

	double a, b;

	uchar *Y_buff_in = new uchar[xSize*ySize];
	char *U_buff_in = new char[xSize / 2 * ySize / 2];
	char *V_buff_in = new char[xSize / 2 * ySize / 2];

	uchar *Y_buff_out = new uchar[xSize*ySize];
	char *U_buff_out = new char[xSize / 2 * ySize / 2];
	char *V_buff_out = new char[xSize / 2 * ySize / 2];

	RGBtoYUV420(input, xSize, ySize, Y_buff_in, U_buff_in, V_buff_in);



	for (int i = 0; i < ySize; i++)
	{
		for (int j = 0; j < xSize; j++)
		{
			jj = j*cos(angle) - i*sin(angle) + m*sin(angle) - n*cos(angle) + n;
			ii = i*cos(angle) + j*sin(angle) - m*cos(angle) - n*sin(angle) + m;				

			if (ii > ySize - 1 || ii < 0 || jj > xSize - 1 || jj < 0)
				Y_buff_out[i*xSize + j] = 0;
			else
			{
				a = jj - floor(jj);
				b = ii - floor(ii);

				m1 = (int) floor(ii);
				m2 = (int) ceil(ii);
				n1 = (int) floor(jj);
				n2 = (int) ceil(jj);

				Y_buff_out[i*xSize + j] = Y_buff_in[m1*xSize + n1] * (1 - a)*(1 - b) +
					Y_buff_in[m2*xSize + n1] * (1 - a)*b
					+ Y_buff_in[m1*xSize + n2] * a*(1 - b)
					+ Y_buff_in[m2*xSize + n2] * a*b;
			}
		}
	}

	for (int i = 0; i < ySize/2; i++)
	{
		for (int j = 0; j < xSize/2; j++)
		{
			jj = j*cos(angle) - i*sin(angle) + m/2.0*sin(angle) - n/2.0*cos(angle) + n/2.0;
			ii = i*cos(angle) + j*sin(angle) - m / 2.0*cos(angle) - n / 2.0*sin(angle) + m / 2.0;

			if (ii > ySize / 2 - 1 || ii < 0 || jj > xSize / 2 - 1 || jj < 0)
			{
				U_buff_out[i*xSize / 2 + j] = 0;
				V_buff_out[i*xSize / 2 + j] = 0;
			}
			else
			{
				a = jj - floor(jj);
				b = ii - floor(ii);

				m1 = (int)floor(ii);
				m2 = (int)ceil(ii);
				n1 = (int)floor(jj);
				n2 = (int)ceil(jj);

				U_buff_out[i*xSize / 2 + j] = U_buff_in[m1*xSize / 2 + n1] * (1 - a)*(1 - b) +
					U_buff_in[m2*xSize / 2 + n1] * (1 - a)*b
					+ U_buff_in[m1*xSize / 2 + n2] * a*(1 - b)
					+ U_buff_in[m2*xSize / 2 + n2] * a*b;

				V_buff_out[i*xSize / 2 + j] = V_buff_in[m1*xSize / 2 + n1] * (1 - a)*(1 - b) +
					V_buff_in[m2*xSize / 2 + n1] * (1 - a)*b
					+ V_buff_in[m1*xSize / 2 + n2] * a*(1 - b)
					+ V_buff_in[m2*xSize / 2 + n2] * a*b;
			}
		}
	}
	YUV420toRGB(Y_buff_out, U_buff_out, V_buff_out, xSize, ySize, output);

	delete Y_buff_in;
	delete U_buff_in;
	delete V_buff_in;
	delete Y_buff_out;
	delete U_buff_out;
	delete V_buff_out;
}