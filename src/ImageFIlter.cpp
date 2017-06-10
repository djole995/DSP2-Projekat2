#include "ImageFilter.h"


void convolve2D (uchar Y_buff[], int xSize, int ySize, double filterCoeff, int N)
{
	//TO DO
}

void extendBorders(uchar input[], int xSize, int ySize, uchar output[], int delta)
{
	for (int i = 0; i < ySize; i++)
	{
		memcpy(&output[i*(xSize + delta)], &input[i*(xSize)], xSize);
		if (delta != 0)
		{
			memset(&output[i*(xSize + delta) + xSize], output[i*(xSize + delta) + xSize - 1], delta);
		}
	}
}
	
void performNFFilter (uchar input[], int xSize, int ySize)
{
	//TO DO
}

void performVFFilter (uchar input[], int xSize, int ySize)
{
	//TO DO
}

void performSuccessiveVFFilter (uchar input[], int xSize, int ySize, int stages)
{
	//TO DO
}

void performSobelEdgeDetection(uchar input[], int xSize, int ySize, uchar threshold)
{
	//TO DO
}

void performNFplusSobelEdgeDetection(uchar input[], int xSize, int ySize, int stages, uchar threshold)
{
	//TO DO
}
