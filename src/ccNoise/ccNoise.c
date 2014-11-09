#include <stdlib.h>

#include <ccNoise/ccNoise.h>
#include <ccRandom/ccRandom.h>
#include <ccTrigonometry/ccTrigonometry.h>

#define CCN_CEIL_INT(n, div) ((n + div - 1) / div)

static float ccnInterpolateLinear(float a, float b, float x)
{
	return a * (1 - x) + b*x;
}

static float ccnInterpolateCosine(float a, float b, float x)
{
	float factor = (1 - ccTriCosDeg(x * 360)) * .5f;

	return a * (1 - factor) + b * factor;
}

static float ccnInterpolateCubic(float a, float b, float c, float d, float x)
{
	float p = (d - c) - (a - b);

	return ccTriCubed(x * p) + ccTriSquared(x * (a - b) - p) + x * (c - a) + b;
}

int ccnGeneratePerlinNoise2D(unsigned int seed, unsigned int width, unsigned int height, unsigned int octaves, unsigned int maxOctave, float persistence, float **buffer)
{
	unsigned int i, j;
	unsigned int bufferSize = width * height;
	unsigned int octaveSize;

	ccrGenerator randomizer;
	unsigned int **randomValues = malloc(sizeof(unsigned int*)* octaves);

	ccrSeed(&randomizer, seed);
	octaveSize = maxOctave;

	for(i = 0; i < octaves; i++) {
		unsigned int octaveWidth = CCN_CEIL_INT(width, octaveSize);
		unsigned int octaveHeight = CCN_CEIL_INT(height, octaveSize);
		unsigned int octaveTotal = octaveWidth * octaveHeight;

		randomValues[i] = malloc(sizeof(unsigned int)* octaveTotal);
		for(j = 0; j < octaveTotal; j++) {
			randomValues[i][j] = ccrGenerateUint(&randomizer);
		}

		octaveSize *= persistence;
	}



	*buffer = malloc(sizeof(float)* bufferSize);




	return 0;
}