#include <ccTrigonometry/ccTrigonometry.h>
#include <stdlib.h>
#include <assert.h>

#include "ccNoiseInternal.h"

#define _CCN_PERLIN_NORMALIZER 0.707107

void ccnGeneratePerlinNoise1D(
	ccnNoise *noise,
	const ccnNoiseConfiguration *configuration,
	const uint32_t scale,
	const ccnInterpolationMethod interpolationMethod)
{
	const uint32_t steps = (uint32_t)ceil((float)noise->width / scale);
	uint32_t interpolationOffset = 0;
	uint32_t period;
	uint32_t i;

	int32_t coordinateOffset = configuration->x * steps;

	float multiplier = configuration->range.high - configuration->range.low;
	float *vectors = malloc(sizeof(float)* (steps + 1));

#ifdef _DEBUG
	assert(interpolationMethod != CCN_INTERP_CUBIC);
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
		period = (uint32_t)(configuration->tileConfiguration.xPeriod * ((float)noise->width / scale));
	}

	for(i = 0; i <= steps; ++i) {
		vectors[i] = cosf((float)(ccrGenerateFloatCoordinate(configuration->seed, ccnWrapCoordinate(i + coordinateOffset, period), 0) * CC_TRI_PI_DOUBLE));
	}

	for(i = 0; i < noise->width; ++i) {
		const uint32_t step = i / scale;
		const float vec = (float)(i + interpolationOffset - step * scale) / scale;

		ccnStore(noise->values + i, configuration->storeMethod, (float)((ccnInterpolate(vectors[step] * vec, vectors[step + 1] * (vec - 1.0f), vec, interpolationMethod) + _CCN_PERLIN_NORMALIZER) * multiplier * _CCN_PERLIN_NORMALIZER) + configuration->range.low);
	}

	free(vectors);
}

void ccnGeneratePerlinNoise2D(
	ccnNoise *noise,
	const ccnNoiseConfiguration *configuration,
	const uint32_t scale,
	const ccnInterpolationMethod interpolationMethod)
{
	const uint32_t xSteps = (uint32_t)ceil((float)noise->width / scale);
	const uint32_t ySteps = (uint32_t)ceil((float)noise->height / scale);
	const uint32_t totalSteps = (xSteps + 1) * (ySteps + 1);
	const uint32_t size = noise->width * noise->height;
	uint32_t xPeriod;
	uint32_t yPeriod;
	uint32_t xOffset = 0;
	uint32_t yOffset = 0;
	uint32_t i;

	float multiplier = configuration->range.high - configuration->range.low;
	float *vectors = malloc(sizeof(float)* (totalSteps << 1));

	ccnPoint offset = (ccnPoint){ configuration->x * xSteps, configuration->y * ySteps };

#ifdef _DEBUG
	assert(interpolationMethod != CCN_INTERP_CUBIC);
	assert(!(scale & (scale - 1)));
#endif

	if(noise->width < scale) {
		offset.x = (int)floor(offset.x * ((float)noise->width / scale));
		xOffset = ccnFloorMod(configuration->x, scale / noise->width) * noise->width;
	}

	if(noise->height < scale) {
		offset.y = (int)floor(offset.y * ((float)noise->height / scale));
		yOffset = ccnFloorMod(configuration->y, scale / noise->height) * noise->height;
	}

	if(configuration->tileConfiguration.tileMethod == CCN_TILE_NOT) {
		xPeriod = yPeriod = CCN_INFINITE;
	}
	else{
		xPeriod = (uint32_t)(configuration->tileConfiguration.xPeriod * ((float)noise->width / scale));
		yPeriod = (uint32_t)(configuration->tileConfiguration.yPeriod * ((float)noise->height / scale));
	}

	for(i = 0; i < totalSteps; ++i) {
		uint32_t Y = i / (xSteps + 1);

		float radians = (float)(ccrGenerateFloatCoordinate(configuration->seed, ccnWrapCoordinate((i - Y * (xSteps + 1)) + offset.x, xPeriod), ccnWrapCoordinate(Y + offset.y, yPeriod)) * CC_TRI_PI_DOUBLE);

		vectors[i << 1] = (float)cos(radians);
		vectors[(i << 1) + 1] = (float)sin(radians);
	}

	for(i = 0; i < size; ++i) {
		uint32_t Y = i / noise->width;
		uint32_t X = i - Y * noise->width;
		uint32_t xStep = X / scale;
		uint32_t yStep = Y / scale;

		X += xOffset;
		Y += yOffset;

		uint32_t indexTop = (xStep + yStep * (xSteps + 1)) << 1;
		uint32_t indexBottom = indexTop + ((xSteps + 1) << 1);
		float vecX = (float)(X - xStep * scale) / scale;
		float vecY = (float)(Y - yStep * scale) / scale;

		ccnStore(noise->values + i, configuration->storeMethod,
			(float)((ccnInterpolate(
			ccnInterpolate(
			vectors[indexTop] * vecX + vectors[indexTop + 1] * vecY,
			vectors[indexTop + 2] * (vecX - 1.0f) + vectors[indexTop + 3] * vecY,
			vecX, interpolationMethod),
			ccnInterpolate(
			vectors[indexBottom] * vecX + vectors[indexBottom + 1] * (vecY - 1.0f),
			vectors[indexBottom + 2] * (vecX - 1.0f) + vectors[indexBottom + 3] * (vecY - 1.0f),
			vecX, interpolationMethod),
			(float)(Y - yStep * scale) / scale, interpolationMethod) + _CCN_PERLIN_NORMALIZER) * _CCN_PERLIN_NORMALIZER * multiplier) + configuration->range.low);
	}

	free(vectors);
}
