#include <stdlib.h>
#include <string.h>

#include <stdio.h>

#include <ccNoise/ccNoise.h>
#include <ccRandom/ccRandom.h>
#include <ccTrigonometry/ccTrigonometry.h>

#define _CCN_CEIL_DIV_INT(n, div) ((n + div - 1) / div)

unsigned int ccnCoordinateUid(int x, int y)
{
	unsigned int shell = max(abs(x), abs(y));
	unsigned int uid = (x + y + (shell << 1)) << 1;
	
	if(shell != 0) {
		uid += ccTriSquared((shell << 1) - 1);
		if(x > y || (x == shell && y == shell)) {
			uid--;
		}
	}

	return uid;
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
}