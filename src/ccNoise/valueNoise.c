#include <math.h>
#include <assert.h>
#include <stdlib.h>

#include "ccNoiseInternal.h"

void ccnGenerateValueNoise1D(
	ccnNoise *noise,
	ccnNoiseConfiguration *configuration,
	unsigned int scale,
	ccnInterpolationMethod interpolationMethod)
{
	unsigned int size = noise->width;
	unsigned int steps = (unsigned int)ceil((float)size / scale);
	unsigned int i, j;
	unsigned int period;
	unsigned int interpolationOffset = 0;
	int coordinateOffset = configuration->x * steps;

	float multiplier = configuration->range.high - configuration->range.low;
	float *bufferedValues;

#ifdef _DEBUG
	assert(!(scale & (scale - 1)));
#endif

	if(noise->width < scale) {
		coordinateOffset = (int)floor(coordinateOffset * ((float)noise->width / scale));
		interpolationOffset = ccnFloorMod(configuration->x, scale / noise->width) * noise->width;
	}

	if(configuration->tileConfiguration.tileMethod == CCN_TILE_NOT) {
		period = CCN_INFINITE;
	}
	else {
		period = (unsigned int)(configuration->tileConfiguration.xPeriod * ((float)noise->width / scale));
	}

	bufferedValues = malloc(sizeof(float)* (interpolationMethod == CCN_INTERP_CUBIC?4:2));

	for(i = 0; i < size; i++) {
		unsigned int oct = i / scale;
		float factor = (float)(i - oct * scale) / scale;

		if(interpolationMethod == CCN_INTERP_CUBIC) {
			if(factor == 0) {
				if(i == 0) {
					for(j = 0; j < 4;) {
						bufferedValues[j] = ccrGenerateFloatCoordinate(configuration->seed, ccnWrapCoordinate(oct - 1 + j++ + coordinateOffset, period), 0);
					}
				}
				else {
					for(j = 0; j < 3;) {
						bufferedValues[j] = bufferedValues[++j];
					}

					bufferedValues[3] = ccrGenerateFloatCoordinate(configuration->seed, ccnWrapCoordinate(oct + 2 + coordinateOffset, period), 0);
				}
			}

			ccnStore(noise->values + i, configuration->storeMethod, (ccnInterpolateCubic(bufferedValues[0], bufferedValues[1], bufferedValues[2], bufferedValues[3], factor + (float)interpolationOffset / scale) * .5f + 0.25f) * multiplier + configuration->range.low);
		}
		else {
			if(factor == 0) {
				if(i == 0) {
					for(j = 0; j < 2;) {
						bufferedValues[j] = ccrGenerateFloatCoordinate(configuration->seed, ccnWrapCoordinate(oct + j++ + coordinateOffset, period), 0);
					}
				}
				else {
					bufferedValues[0] = bufferedValues[1];
				}

				bufferedValues[1] = ccrGenerateFloatCoordinate(configuration->seed, ccnWrapCoordinate(oct + 1 + coordinateOffset, period), 0);
			}

			ccnStore(noise->values + i, configuration->storeMethod, ccnInterpolate(bufferedValues[0], bufferedValues[1], factor + (float)interpolationOffset / scale, interpolationMethod) * multiplier + configuration->range.low);
		}
	}

	free(bufferedValues);
}

void ccnGenerateValueNoise2D(
	ccnNoise *noise,
	ccnNoiseConfiguration *configuration,
	unsigned int scale,
	ccnInterpolationMethod interpolationMethod)
{
	unsigned int size = noise->width * noise->height;
	unsigned int xSteps = (unsigned int)ceil((float)noise->width / scale);
	unsigned int ySteps = (unsigned int)ceil((float)noise->height / scale);
	unsigned int offsetHeight;
	unsigned int xPeriod;
	unsigned int yPeriod;
	unsigned int xOffset = 0;
	unsigned int yOffset = 0;
	unsigned int interpolationYoffset = interpolationMethod == CCN_INTERP_CUBIC?1:0;
	unsigned int i, j, k;

	float multiplier = configuration->range.high - configuration->range.low;
	float *xValues;
	float *bufferedValues;

	ccnPoint coordinateOffset = (ccnPoint){ configuration->x * xSteps, configuration->y * ySteps };

#ifdef _DEBUG
	assert(!(scale & (scale - 1)));
#endif

	if(noise->width < scale) {
		coordinateOffset.x = (int)floor(coordinateOffset.x * ((float)noise->width / scale));
		xOffset = ccnFloorMod(configuration->x, scale / noise->width) * noise->width;
	}

	if(noise->height < scale) {
		coordinateOffset.y = (int)floor(coordinateOffset.y * ((float)noise->height / scale));
		yOffset = ccnFloorMod(configuration->y, scale / noise->height) * noise->height;
	}

	if(configuration->tileConfiguration.tileMethod == CCN_TILE_NOT) {
		xPeriod = yPeriod = CCN_INFINITE;
	}
	else{
		xPeriod = (unsigned int)(configuration->tileConfiguration.xPeriod * ((float)noise->width / scale));
		yPeriod = (unsigned int)(configuration->tileConfiguration.yPeriod * ((float)noise->height / scale));
	}

	offsetHeight = ySteps + (interpolationMethod == CCN_INTERP_CUBIC?3:1);
	xValues = malloc(noise->width * offsetHeight * sizeof(float));
	bufferedValues = malloc(sizeof(float)* (interpolationMethod == CCN_INTERP_CUBIC?4:2));

	for(i = 0; i < offsetHeight; ++i) {
		for(j = 0; j < noise->width; ++j) {
			unsigned int octX = j / scale;

			float factor = (float)(j - octX * scale) / scale;

			if(interpolationMethod == CCN_INTERP_CUBIC) {
				if(factor == 0) {
					if(j == 0) {
						for(k = 0; k < 4; k++) {
							bufferedValues[k] = ccrGenerateFloatCoordinate(configuration->seed, ccnWrapCoordinate(octX - 1 + k + coordinateOffset.x, xPeriod), ccnWrapCoordinate(i + coordinateOffset.y, yPeriod));
						}
					}
					else {
						for(k = 0; k < 3; k++) {
							bufferedValues[k] = bufferedValues[k + 1];
						}

						bufferedValues[3] = ccrGenerateFloatCoordinate(configuration->seed, ccnWrapCoordinate(octX + 2 + coordinateOffset.x, xPeriod), ccnWrapCoordinate(i + coordinateOffset.y, yPeriod));
					}
				}

				xValues[i * noise->width + j] = ccnInterpolateCubic(bufferedValues[0], bufferedValues[1], bufferedValues[2], bufferedValues[3], factor + (float)xOffset / scale);
			}
			else {
				if(factor == 0) {
					if(j == 0) {
						bufferedValues[0] = ccrGenerateFloatCoordinate(configuration->seed, ccnWrapCoordinate(octX + coordinateOffset.x, xPeriod), ccnWrapCoordinate(i + coordinateOffset.y, yPeriod));
					}
					else {
						bufferedValues[0] = bufferedValues[1];
					}

					bufferedValues[1] = ccrGenerateFloatCoordinate(configuration->seed, ccnWrapCoordinate(octX + 1 + coordinateOffset.x, xPeriod), ccnWrapCoordinate(i + coordinateOffset.y, yPeriod));
				}

				xValues[i * noise->width + j] = ccnInterpolate(bufferedValues[0], bufferedValues[1], factor + (float)xOffset / scale, interpolationMethod);
			}
		}
	}

	for(i = 0; i < size; ++i) {
		unsigned int Y = i / noise->width;
		unsigned int octY = Y / scale;
		unsigned int index = (i - Y * noise->width) + (octY + interpolationYoffset) * noise->width;

		float factor = (float)(Y - octY * scale) / scale;
		float interpFactor = factor + (float)yOffset / scale;

		if(interpolationMethod == CCN_INTERP_CUBIC) {
			ccnStore(noise->values + i, configuration->storeMethod, (ccnInterpolateCubic(xValues[index - noise->width], xValues[index], xValues[index + noise->width], xValues[index + (noise->width << 1)], interpFactor) * .5f + 0.25f) * multiplier + configuration->range.low);
		}
		else {
			ccnStore(noise->values + i, configuration->storeMethod, ccnInterpolate(xValues[index], xValues[index + noise->width], interpFactor, interpolationMethod) * multiplier + configuration->range.low);
		}
	}

	free(xValues);
	free(bufferedValues);
}