#include <stdlib.h>

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

int ccnGeneratePerlinNoise2D(unsigned int seed, unsigned int width, unsigned int height, unsigned int octaves, unsigned int maxOctave, float octavePersistence, float **buffer2)
{
	unsigned int i, j;
	unsigned int bufferSize = width * height;
	unsigned int octaveSize;
	float factor = 0.5f;

	ccrGenerator randomizer;

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
		factor /= 2;
	}

	*buffer2 = buffer;

	return 0;
}