#include <stdlib.h>
#include <string.h>

#include <stdio.h>
#include <stdlib.h>

#include <ccNoise/ccNoise.h>
#include <ccRandom/ccRandom.h>
#include <ccTrigonometry/ccTrigonometry.h>
#include <ccSort/ccSort.h>

static float ccnInterpolate(float a, float b, float x, ccnInterpolationMethod interpolationMethod)
{
	switch(interpolationMethod) {
	case CCN_INTERP_LINEAR:
		return ccTriInterpolateLinear(a, b, x);
	case CCN_INTERP_QUADRATIC:
		return ccTriInterpolateQuadratic(a, b, x);
	case CCN_INTERP_QUADRATIC_INVERSE:
		return ccTriInterpolateQuadraticInverse(a, b, x);
	case CCN_INTERP_COSINE:
		return ccTriInterpolateCosine(a, b, x);
	default:
		return 0;
	}
}

static unsigned int ccnDistance(ccPoint a, ccPoint b, ccnDistanceMethod distanceMethod)
{
	switch(distanceMethod) {
	case CCN_DIST_MANHATTAN:
		return abs(a.x - b.x) + abs(a.y - b.y);
	case CCN_DIST_EUCLIDEAN:
		return (unsigned int)ccTriDistance(a.x, a.y, b.x, b.y);
	case CCN_DIST_CHEBYCHEV:
		return max(abs(a.x - b.x), abs(a.y - b.y));
	default:
		return 0;
	}
}

static unsigned int ccnGenerateNoiseSeed(unsigned int seed, int x, int y)
{
	return seed + ccnCoordinateUid(x, y);
}

static int ccnWrapCoordinate(int coordinate, unsigned int period) {
	int positiveDeviation, negativeDeviation;

	if(period == CCN_INFINITE) return coordinate;

	positiveDeviation = period >> 1;
	negativeDeviation = -(int)(period - positiveDeviation);

	if(coordinate > positiveDeviation) {
		coordinate -= ((coordinate - negativeDeviation - 1) / period) * period;
	}
	else if(coordinate <= negativeDeviation) {
		coordinate += ((-coordinate + positiveDeviation) / period) * period;
	}

	return coordinate;
}

static void ccnGenerateOffsetNoise(
	float **buffer,
	unsigned int seed,
	ccnTileConfiguration *tileConfig,
	int x, int y,
	unsigned int width, unsigned int height,
	ccPoint negativeOffset, ccPoint positiveOffset)
{
	unsigned int X, Y;

	unsigned int totalWidth = width + negativeOffset.x + positiveOffset.x;
	unsigned int totalHeight = height + negativeOffset.y + positiveOffset.y;

	if(tileConfig->tileMethod == CCN_TILE_NOT) {
		ccnGenerateWhiteNoise(buffer, ccnGenerateNoiseSeed(seed, x, y), totalWidth, totalHeight);
	}
	else {
		float *whiteNoiseBuffer = malloc(sizeof(float) * totalWidth * totalHeight);
		bool wrapCorners = tileConfig->xPeriod == 1 && tileConfig->yPeriod == 1;

		*buffer = malloc(sizeof(float)*totalWidth * totalHeight);

		// Generate central noise
		ccnGenerateWhiteNoise(&whiteNoiseBuffer, ccnGenerateNoiseSeed(seed, ccnWrapCoordinate(x, tileConfig->xPeriod), ccnWrapCoordinate(y, tileConfig->yPeriod)), width, height);

		for(Y = 0; Y < height; Y++) {
			for(X = 0; X < width; X++) {
				(*buffer)[(Y + negativeOffset.y) * totalWidth + X + negativeOffset.x] = whiteNoiseBuffer[(Y * width) + X];
			}
		}
		
		// Generate horizontal noises

		// Positive
		if(positiveOffset.x > 0) {

#define _CCN_BUFFERDEST (*buffer)[-positiveOffset.x + X + totalWidth * (negativeOffset.y + Y + 1)]

			if(tileConfig->xPeriod == 1) {
				for(Y = 0; Y < height; Y++) {
					for(X = 0; X < (unsigned int)positiveOffset.x; X++) {
						_CCN_BUFFERDEST = (*buffer)[negativeOffset.x + X + totalWidth * (negativeOffset.y + Y)];
					}
				}
			}
			else {
				ccnGenerateWhiteNoise(&whiteNoiseBuffer, ccnGenerateNoiseSeed(seed, ccnWrapCoordinate(x + 1, tileConfig->xPeriod), ccnWrapCoordinate(y, tileConfig->yPeriod)), width, height);

				for(Y = 0; Y < height; Y++) {
					for(X = 0; X < (unsigned int)positiveOffset.x; X++) {
						_CCN_BUFFERDEST = whiteNoiseBuffer[X + Y * width];
					}
				}
			}
		}

		// Negative
		if(negativeOffset.x > 0) {

#undef _CCN_BUFFERDEST
#define _CCN_BUFFERDEST (*buffer)[X + totalWidth * (negativeOffset.y + Y)]

			if(tileConfig->xPeriod == 1) {
				for(Y = 0; Y < height; Y++) {
					for(X = 0; X < (unsigned int)negativeOffset.x; X++) {
						_CCN_BUFFERDEST = (*buffer)[width + X + totalWidth * (Y + negativeOffset.y)];
					}
				}
			}
			else {
				ccnGenerateWhiteNoise(&whiteNoiseBuffer, ccnGenerateNoiseSeed(seed, ccnWrapCoordinate(x - 1, tileConfig->xPeriod), ccnWrapCoordinate(y, tileConfig->yPeriod)), width, height);

				for(Y = 0; Y < height; Y++) {
					for(X = 0; X < (unsigned int)negativeOffset.x; X++) {
						_CCN_BUFFERDEST = whiteNoiseBuffer[-negativeOffset.x + X + width * (Y + 1)];
					}
				}
			}
		}

		// Generate vertical noises

		// Positive
		if(positiveOffset.y > 0) {

#undef _CCN_BUFFERDEST
#define _CCN_BUFFERDEST (*buffer)[negativeOffset.x + X + totalWidth * (height + negativeOffset.y + Y)]

			if(tileConfig->yPeriod == 1) {
				for(X = 0; X < width; X++) {
					for(Y = 0; Y < (unsigned int)positiveOffset.y; Y++) {
						_CCN_BUFFERDEST = (*buffer)[negativeOffset.x + X + totalWidth * (negativeOffset.y + Y)];
					}
				}
			}
			else {
				ccnGenerateWhiteNoise(&whiteNoiseBuffer, ccnGenerateNoiseSeed(seed, ccnWrapCoordinate(x, tileConfig->xPeriod), ccnWrapCoordinate(y + 1, tileConfig->yPeriod)), width, positiveOffset.y);

				for(X = 0; X < width; X++) {
					for(Y = 0; Y < (unsigned int)positiveOffset.y; Y++) {
						_CCN_BUFFERDEST = whiteNoiseBuffer[X + width * Y];
					}
				}
			}
		}

		// Negative
		if(negativeOffset.y > 0) {

#undef _CCN_BUFFERDEST
#define _CCN_BUFFERDEST (*buffer)[negativeOffset.x + X + totalWidth * Y]

			if(tileConfig->yPeriod == 1) {
				for(X = 0; X < width; X++) {
					for(Y = 0; Y < (unsigned int)negativeOffset.y; Y++) {
						_CCN_BUFFERDEST = (*buffer)[negativeOffset.x + X + totalWidth * (height + Y)];
					}
				}
			}
			else {
				ccnGenerateWhiteNoise(&whiteNoiseBuffer, ccnGenerateNoiseSeed(seed, ccnWrapCoordinate(x, tileConfig->xPeriod), ccnWrapCoordinate(y - 1, tileConfig->yPeriod)), width, height);

				for(X = 0; X < width; X++) {
					for(Y = 0; Y < (unsigned int)negativeOffset.y; Y++) {
						_CCN_BUFFERDEST = whiteNoiseBuffer[X + width * (height - negativeOffset.y + Y)];
					}
				}
			}
		}

		// Generate corner noises

		// Right bottom
		if(positiveOffset.x > 0 && positiveOffset.y > 0) {

#undef _CCN_BUFFERDEST
#define _CCN_BUFFERDEST (*buffer)[width + negativeOffset.x + X + totalWidth * (height + negativeOffset.y + Y)]

			if(wrapCorners) {
				for(X = 0; X < (unsigned int)positiveOffset.x; X++) {
					for(Y = 0; Y < (unsigned int)positiveOffset.y; Y++) {
						_CCN_BUFFERDEST = (*buffer)[negativeOffset.x + X + totalWidth * (negativeOffset.y + Y)];
					}
				}
			}
			else {
				ccnGenerateWhiteNoise(&whiteNoiseBuffer, ccnGenerateNoiseSeed(seed, ccnWrapCoordinate(x + 1, tileConfig->xPeriod), ccnWrapCoordinate(y + 1, tileConfig->yPeriod)), width, positiveOffset.y);

				for(X = 0; X < (unsigned int)positiveOffset.x; X++) {
					for(Y = 0; Y < (unsigned int)positiveOffset.y; Y++) {
						_CCN_BUFFERDEST = whiteNoiseBuffer[X + width * Y];
					}
				}
			}
		}

		// Left bottom
		if(negativeOffset.x > 0 && positiveOffset.y > 0) {

#undef _CCN_BUFFERDEST
#define _CCN_BUFFERDEST (*buffer)[X + totalWidth * (height + negativeOffset.y + Y)]

			if(wrapCorners) {
				for(X = 0; X < (unsigned int)negativeOffset.x; X++) {
					for(Y = 0; Y < (unsigned int)positiveOffset.y; Y++) {
						_CCN_BUFFERDEST = (*buffer)[width + X + totalWidth * (negativeOffset.y + Y)];
					}
				}
			}
			else {
				ccnGenerateWhiteNoise(&whiteNoiseBuffer, ccnGenerateNoiseSeed(seed, ccnWrapCoordinate(x - 1, tileConfig->xPeriod), ccnWrapCoordinate(y + 1, tileConfig->yPeriod)), width, positiveOffset.y);

				for(X = 0; X < (unsigned int)negativeOffset.x; X++) {
					for(Y = 0; Y < (unsigned int)positiveOffset.y; Y++) {
						_CCN_BUFFERDEST = whiteNoiseBuffer[-negativeOffset.x + X + width * (Y + 1)];
					}
				}
			}
		}

		// Left top
		if(negativeOffset.x > 0 && negativeOffset.y > 0) {

#undef _CCN_BUFFERDEST
#define _CCN_BUFFERDEST (*buffer)[X + totalWidth * Y]

			if(wrapCorners) {
				for(X = 0; X < (unsigned int)negativeOffset.x; X++) {
					for(Y = 0; Y < (unsigned int)negativeOffset.y; Y++) {
						_CCN_BUFFERDEST = (*buffer)[width + X + totalWidth * (height + Y)];
					}
				}
			}
			else {
				ccnGenerateWhiteNoise(&whiteNoiseBuffer, ccnGenerateNoiseSeed(seed, ccnWrapCoordinate(x - 1, tileConfig->xPeriod), ccnWrapCoordinate(y - 1, tileConfig->yPeriod)), width, height);

				for(X = 0; X < (unsigned int)negativeOffset.x; X++) {
					for(Y = 0; Y < (unsigned int)negativeOffset.y; Y++) {
						_CCN_BUFFERDEST = whiteNoiseBuffer[-negativeOffset.x + X + width * (height - negativeOffset.y + Y + 1)];
					}
				}
			}
		}

		// Right top
		if(positiveOffset.x > 0 && negativeOffset.y > 0) {

#undef _CCN_BUFFERDEST
#define _CCN_BUFFERDEST (*buffer)[width + negativeOffset.x + X + totalWidth * Y]

			if(wrapCorners) {
				for(X = 0; X < (unsigned int)positiveOffset.x; X++) {
					for(Y = 0; Y < (unsigned int)negativeOffset.y; Y++) {
						_CCN_BUFFERDEST = (*buffer)[negativeOffset.x + X + totalWidth * (height + Y)];
					}
				}
			}
			else {
				ccnGenerateWhiteNoise(&whiteNoiseBuffer, ccnGenerateNoiseSeed(seed, ccnWrapCoordinate(x + 1, tileConfig->xPeriod), ccnWrapCoordinate(y - 1, tileConfig->yPeriod)), width, height);

				for(X = 0; X < (unsigned int)positiveOffset.x; X++) {
					for(Y = 0; Y < (unsigned int)negativeOffset.y; Y++) {
						_CCN_BUFFERDEST = whiteNoiseBuffer[X + width * (height - negativeOffset.y + Y)];
					}
				}
			}
		}

		free(whiteNoiseBuffer);
	}
}

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

int ccnGenerateWorleyNoise(
	float **buffer,
	unsigned int seed,
	ccnTileConfiguration *tileConfig,
	int x, int y,
	unsigned int width, unsigned int height,
	unsigned int points,
	unsigned int n,
	int low, int high,
	float lowValue, float highValue,
	ccnDistanceMethod distanceMethod,
	ccnInterpolationMethod interpolationMethod)
{
	unsigned int size = width * height;
	unsigned int i, j;
	unsigned int pointId = 0;
	unsigned int pointListSize = points * 9;
	ccPoint offset;

	ccPoint *pointList = malloc(pointListSize*sizeof(ccPoint));
	int *pointsDistances = malloc(pointListSize*sizeof(unsigned int));

	unsigned int maxManhattanDistance = (unsigned int)(high * (2 / sqrt(2)));

	if(interpolationMethod == CCN_INTERP_CUBIC) return CCN_ERROR_INVALID_METHOD;

	if(tileConfig->tileMethod = CCN_TILE_NOT) tileConfig->xPeriod = tileConfig->yPeriod = CCN_INFINITE;

	for(offset.x = -1; offset.x <= 1; offset.x++) {
		for(offset.y = -1; offset.y <= 1; offset.y++) {
			ccRandomizer32 randomizer;

			ccrSeed32(&randomizer, ccnGenerateNoiseSeed(seed, ccnWrapCoordinate(x + offset.x, tileConfig->xPeriod), ccnWrapCoordinate(y + offset.y, tileConfig->yPeriod)));

			for(i = 0; i < points; i++) {
				pointList[pointId].x = (int)((ccrGenerateFloat32(&randomizer) + offset.x) * width);
				pointList[pointId].y = (int)((ccrGenerateFloat32(&randomizer) + offset.y) * height);
				pointId++;
			}
		}
	}

	for(i = 0; i < size; i++) {
		ccPoint p;
		p.y = i / width;
		p.x = i - p.y * width;

		pointId = 0;
		for(j = 0; j < pointListSize; j++) {
			unsigned int manhattanDistance = ccnDistance(p, pointList[j], CCN_DIST_MANHATTAN);

			if(manhattanDistance < maxManhattanDistance) {
				if(distanceMethod == CCN_DIST_MANHATTAN) {
					pointsDistances[pointId] = manhattanDistance;
				}
				else {
					pointsDistances[pointId] = ccnDistance(p, pointList[j], distanceMethod);
				}

				pointId++;
			}
		}

		if(pointId > 1) ccsQuicksort(pointsDistances, 0, pointId);
		
		if(pointId <= n) {
			(*buffer)[i] = highValue;
		}
		else if(pointsDistances[n] > high) {
			(*buffer)[i] = highValue;
		}
		else if(pointsDistances[n] < low) {
			(*buffer)[i] = lowValue;
		}
		else {
			(*buffer)[i] = ccnInterpolate(lowValue, highValue, (float)(pointsDistances[n] - low) / (high - low), interpolationMethod);
		}
	}

	free(pointList);
	free(pointsDistances);

	return CCN_ERROR_NONE;
}

int ccnGenerateWhiteNoise(
	float **buffer,
	unsigned int seed,
	unsigned int width, unsigned int height)
{
	unsigned int size = width * height;
	unsigned int i;
	
	ccRandomizer32 randomizer;

	ccrSeed32(&randomizer, seed);

	for(i = 0; i < size; i++) {
		(*buffer)[i] = ccrGenerateFloat32(&randomizer);
	}

	return CCN_ERROR_NONE;
}

int ccnGenerateValueNoise(
	float **buffer,
	unsigned int seed,
	ccnTileConfiguration *tileConfig,
	int x, int y,
	unsigned int width, unsigned int height,
	unsigned int octaves,
	unsigned int maxOctave,
	ccnInterpolationMethod interpolationMethod)
{
	unsigned int size = width * height;
	unsigned int octaveSize = maxOctave;
	unsigned int i, j, k;
	
	float influence = 0.5f;

	ccPoint negativeOffset;
	ccPoint positiveOffset;

	if(interpolationMethod == CCN_INTERP_CUBIC) {
		negativeOffset.x = negativeOffset.y = 1;
		positiveOffset.x = positiveOffset.y = 2;
	}
	else {
		negativeOffset.x = negativeOffset.y = 0;
		positiveOffset.x = positiveOffset.y = 1;
	}

	octaveSize = maxOctave;

	for(i = 0; i < size; i++) {
		(*buffer)[i] = 0;
	}

	for(i = 0; i < octaves; i++) {
		unsigned int xSteps = width / octaveSize;
		unsigned int ySteps = height / octaveSize;
		unsigned int offsetWidth = xSteps + negativeOffset.x + positiveOffset.x;
		unsigned int offsetHeight = ySteps + negativeOffset.y + positiveOffset.y;

		float *offsetNoise;
		float *xValues = malloc(width * offsetHeight * sizeof(float));

		ccnGenerateOffsetNoise(&offsetNoise, seed, tileConfig, x, y, xSteps, ySteps, negativeOffset, positiveOffset);

		for(j = 0; j < offsetHeight; j++) {
			for(k = 0; k < width; k++) {
				unsigned int octX = k / octaveSize;
				unsigned int offsetIndex = j * offsetWidth + negativeOffset.x + octX;

				float factor = (float)(k - octX * octaveSize) / octaveSize;

				if(factor == 0) {
					xValues[j * width + k] = offsetNoise[offsetIndex];
				}
				else {
					if(interpolationMethod == CCN_INTERP_CUBIC) {
						xValues[j * width + k] = ccTriInterpolateCubic(offsetNoise[offsetIndex - 1], offsetNoise[offsetIndex], offsetNoise[offsetIndex + 1], offsetNoise[offsetIndex + 2], factor);
					}
					else {
						xValues[j * width + k] = ccnInterpolate(offsetNoise[offsetIndex], offsetNoise[offsetIndex + 1], factor, interpolationMethod);
					}
				}
			}
		}

		for(j = 0; j < size; j++) {
			unsigned int Y = j / width;
			unsigned int X = j - Y * width;

			unsigned octY = Y / octaveSize;

			float factor = (float)(Y - octY * octaveSize) / octaveSize;

			unsigned int index = X + (octY + negativeOffset.y) * width;

			if(factor == 0) {
				(*buffer)[j] += xValues[index] * influence;
			}
			else {
				if(interpolationMethod == CCN_INTERP_CUBIC) {
					(*buffer)[j] += ccTriInterpolateCubic(xValues[index - width], xValues[index], xValues[index + width], xValues[index + (width << 1)], factor) * influence;
				}
				else {
					(*buffer)[j] += ccnInterpolate(xValues[index], xValues[index + width], factor, interpolationMethod) * influence;
				}
			}
		}
		
		free(xValues);
		free(offsetNoise);

		influence *= 0.5f;

		octaveSize >>= 1;
		if(octaveSize == 0) {
			if(octaves == CCN_INFINITE) {
				break;
			}
			else {
				return CCN_ERROR_INVALID_ARGUMENT_RANGE;
			}
		}
	}

	return CCN_ERROR_NONE;
}