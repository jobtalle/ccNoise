#include <stdlib.h>
#include <string.h>

#include <stdio.h>

#include <ccNoise/ccNoise.h>
#include <ccRandom/ccRandom.h>
#include <ccTrigonometry/ccTrigonometry.h>

#define _CCN_CEIL_DIV_INT(n, div) ((n + div - 1) / div)

typedef struct {
	int x, y;
} ccnPoint;

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

void ccnGenerateWorleyNoise(float **buffer, unsigned int seed, int x, int y, unsigned int width, unsigned int height, unsigned int points)
{
	unsigned int size = width * height;
	unsigned int surface;
	unsigned int i;
	unsigned int pointId = 0;
	ccnPoint *pointList = malloc((points << 2)*sizeof(ccnPoint));

	*buffer = malloc(sizeof(float)*size);

	for(surface = 0; surface < 4; surface++) {
		ccRandomizer32 randomizer;
		ccnPoint offset;

		switch(surface) {
		case 0:
			offset = (ccnPoint){ 0, 0 };
			break;
		case 1:
			offset = (ccnPoint){ 1, 0 };
			break;
		case 2:
			offset = (ccnPoint){ 0, 1 };
			break;
		case 3:
			offset = (ccnPoint){ 1, 1 };
			break;
		}

		ccrSeed32(&randomizer, seed + ccnCoordinateUid(x + offset.x, y + offset.y));

		for(i = 0; i < points; i++) {
			pointList[pointId].x = (unsigned int)((ccrGenerateFloat32(&randomizer) + offset.x) * width);
			pointList[pointId].y = (unsigned int)((ccrGenerateFloat32(&randomizer) + offset.y) * height);
			pointId++;
		}
	}

	for(i = 0; i < size; i++) {
		(*buffer)[i] = 0.5f;
	}
}