#include <stdlib.h>
#include <string.h>

#include <stdio.h>

#include <ccNoise/ccNoise.h>
#include <ccRandom/ccRandom.h>
#include <ccTrigonometry/ccTrigonometry.h>

#define _CCN_CEIL_DIV_INT(n, div) ((n + div - 1) / div)

uint64_t ccnCoordinateUid(int x, int y)
{
	unsigned int shell = max(abs(x), abs(y));
	uint64_t uid = (x + y + (shell << 1)) << 1;
	
	if(shell != 0) {
		uid += ccTriSquared((shell << 1) - 1);
		if(x > y || (x == shell && y == shell)) {
			uid--;
		}
	}

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