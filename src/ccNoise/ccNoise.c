#include <stdlib.h>
#include <string.h>

#include <stdio.h>

#include <ccNoise/ccNoise.h>
#include <ccRandom/ccRandom.h>
#include <ccTrigonometry/ccTrigonometry.h>

#define _CCN_CEIL_DIV_INT(n, div) ((n + div - 1) / div)

static float _ccnInterpolateLinear(float a, float b, float x)
{
	return a * (1 - x) + b*x;
}

static float _ccnInterpolateCosine(float a, float b, float x)
{
	float factor = (1.0f - ccTriCosDeg((unsigned int)(x * 180.0f))) * .5f;

	return a * (1.0f - factor) + b * factor;
}

static float _ccnInterpolateCubic(float a, float b, float c, float d, float x)
{
	float p = (d - c) - (a - b);

	return ccTriCubed(x * p) + ccTriSquared(x * (a - b) - p) + x * (c - a) + b;
}

static unsigned int _ccnAbsMax(unsigned int elementCount, unsigned int *elements)
{
	unsigned int i;
	unsigned int max = abs(elements[0]);

	for(i = 1; i < elementCount; i++) {
		if((unsigned int)abs(elements[i]) > max) {
			max = abs(elements[i]);
		}
	}

	return max;
}

static unsigned int _ccnIntPow(unsigned int base, unsigned int power)
{
	unsigned int i;
	unsigned int result = 0;
	
	for(i = 0; i < power; i++) {
		if(result == 0) {
			result = base;
		}
		else {
			result *= base;
		}
	}

	return result;
}

unsigned int ccnCoordinateUid(unsigned int dimensions, int *coordinate)
{
	int i;
	unsigned int uid;
	unsigned int shell = _ccnAbsMax(dimensions, coordinate);
	unsigned int shellDiameter = ((shell - 1) << 1) + 1;
	
	if(shell == 0) return 0;

	uid = _ccnIntPow(shellDiameter, dimensions);
	/*
	for(i = dimensions - 1; i >= 0; i--) {
		uid += (coordinate[i] + shell) * (_ccnIntPow(shellDiameter, i) - _ccnIntPow(shellDiameter - 2, i));
	}
	*/
	//uid += coordinate[0] + shell;

	return uid;
}

static void _ccnGenerateValueNoise(ccRandomizer *randomizer, float *buffer, unsigned int dimension, unsigned int *sizes, unsigned int octave, unsigned int octaveSize)
{
	/*
	float influence;
	unsigned int i, j;
	unsigned int octaveSize;

	if(dimension == 1) {
		// Random function

		// Return interpolated 1D perlin noise

		octaveSize = maxOctave;
		influence = 0.5f; // TODO: make dynamic

		for(i = 0; i < octaves; i++) {
			unsigned int octaveIntervals = _CCN_CEIL_DIV_INT(sizes[0], octaveSize) + 1;
			float *randomValues = malloc(sizeof(float)* octaveIntervals);

			octaveSize = (unsigned int)((float)octaveSize * octavePersistence);

			for(j = 0; j < octaveIntervals; j++) {
				randomValues[j] = ccrGenerateFloat(randomizer);
			}

			for(j = 0; j < sizes[0]; j++) {
				unsigned int left, right;

				left = j / octaveSize;
				right = left + 1;

				buffer[j] += influence * _ccnInterpolateCosine(randomValues[left], randomValues[right], (float)(j % octaveSize) / octaveSize);
			}

			influence *= octavePersistence;
			(float)octaveSize *= octavePersistence;

			if(octaveSize == 0) {
				printf("Octave size reached zero!\n");
				break;
			}

			free(randomValues);
		}
	}
	else {
		// Recursive calls
		
		// Return vertical interpolations between YFIT 1D noises

		for(i = 0; i < octaves; i++) {
			printf("Octave\n");
		}
	}
	*/
}

// TODO: octave -1 means iterate until smallest detail
void ccnGenerateValueNoise(unsigned int seed, float **buffer, unsigned int dimensions, unsigned int *sizes, unsigned int octaves, unsigned int maxOctave, float octavePersistence)
{
	unsigned int bufferSize = 1;
	unsigned int i, j;
	unsigned int octaveSize = maxOctave;
	float *tempBuffer;
	float influence = 0.5f;

	for(i = 0; i < dimensions; i++) {
		bufferSize *= sizes[i];
	}

	*buffer = calloc(bufferSize, sizeof(float));
	tempBuffer = malloc(bufferSize * sizeof(float));

	ccRandomizer randomizer;
	ccrSeed(&randomizer, seed);

	for(i = 0; i < octaves; i++) {
		printf("Octave %d...\n", i);

		memset(tempBuffer, 0, bufferSize * sizeof(float));

		_ccnGenerateValueNoise(&randomizer, tempBuffer, dimensions, sizes, i, octaveSize);

		for(j = 0; j < bufferSize; j++) {
			// Add temp to buffer with multiplication
			(*buffer)[j] += tempBuffer[j] * influence;
		}

		if(octaveSize == 0 || i == octaves - 1) {
			// Last octave reached
			break;
		}
		else {
			// Prepare for next round
			(float)octaveSize *= octavePersistence;
			influence *= .5f;
		}
	}

	free(tempBuffer);
}

/* Ye olde function
float *ccnGeneratePerlinNoise(unsigned int seed, unsigned int dimensions, unsigned int *sizes, unsigned int octaves, unsigned int maxOctave, float octavePersistence)
{
	unsigned int i, j;
	unsigned int bufferSize = width * height;
	unsigned int octaveSize;
	float factor = 0.5f;

	ccRandomizer randomizer;

	float *buffer = calloc(bufferSize, sizeof(float));
	
	ccrSeed(&randomizer, seed);
	octaveSize = maxOctave;
	
	for(i = 0; i < octaves; i++) { // TODO: octaves
		unsigned int octavesWidth = _CCN_CEIL_DIV_INT(width, octaveSize) + 1;
		unsigned int octavesHeight = _CCN_CEIL_DIV_INT(height, octaveSize) + 1;
		unsigned int octavesTotal = octavesWidth * octavesHeight;

		unsigned int x, y;
		float *randomValues = malloc(sizeof(float)* octavesTotal);

		for(j = 0; j < octavesTotal; j++) {
			randomValues[j] = (float)(ccrGenerateUint(&randomizer) / 4294967295.0f); // TODO: implement ccrGenerateFloat
		}

		x = y = 0;
		for(j = 0; j < bufferSize; j++) {
			unsigned int left, right;
			float topVal, bottomVal;

			left = x / octaveSize;
			right = left + 1;
			
			topVal = _ccnInterpolateCosine(randomValues[left + (y / octaveSize) * octavesWidth], randomValues[right + (y / octaveSize) * octavesWidth], (float)((float)(x % octaveSize) / octaveSize));
			bottomVal = _ccnInterpolateCosine(randomValues[left + ((y / octaveSize) + 1) * octavesWidth], randomValues[right + ((y / octaveSize) + 1) * octavesWidth], (float)((float)(x % octaveSize) / octaveSize));

			buffer[j] += _ccnInterpolateCosine(topVal, bottomVal, (float)((float)(y % octaveSize) / octaveSize)) * factor;

			x++;
			if(x == width) {
				x = 0;
				y++;
			}
		}

		free(randomValues);

		octaveSize = (unsigned int)((float)octaveSize * octavePersistence);
		if(octaveSize == 0) break;
		factor *= 0.5f;
	}

	return buffer;
}
*/