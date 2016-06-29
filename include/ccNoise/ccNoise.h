//__________________________________________________________________________________//
//                                 _   _       _                                    //
//                                | \ | |     (_)                                   //
//                         ___ ___|  \| | ___  _ ___  ___                           //
//                        / __/ __| . ` |/ _ \| / __|/ _ \                          //
//                       | (_| (__| |\  | (_) | \__ \  __/                          //
//                        \___\___|_| \_|\___/|_|___/\___| 1.0                      //
//                                                                                  //
//              Copyright (C) 2014 - 2015 \ Job Talle (jobtalle@hotmail.com)        //
//__________________________________________________________________________________//
//                                                                                  //
//      This program is free software: you can redistribute it and/or modify        //
//      it under the terms of the 3-clause BSD license.                             //
//                                                                                  //
//      You should have received a copy of the 3-clause BSD License along with      //
//      this program. If not, see <http://opensource.org/licenses/>.                //
//__________________________________________________________________________________//

#pragma once

#include <stdint.h>

#ifdef __cplusplus
extern "C"
{
#endif

#define CCN_INFINITE UINT32_MAX

#define ccnNoiseAllocate1D(noise, w) noise.values = malloc(sizeof(float)* w); noise.width = w;
#define ccnNoiseAllocate2D(noise, w, h) noise.values = malloc(sizeof(float)* w * h); noise.width = w; noise.height = h;
#define ccnNoiseGet1D(noise, x) noise.values[x]
#define ccnNoiseGet2D(noise, x, y) noise.values[x + noise.line * y]
#define ccnNoiseSet1D(noise, x, value) noise.values[x] = value;
#define ccnNoiseSet2D(noise, x, y, value) noise.values[x + noise.line * y] = value
#define ccnNoiseFree(noise) free(noise.values)

typedef enum {
	CCN_INTERP_LINEAR,
	CCN_INTERP_COSINE,
	CCN_INTERP_CUBIC,
	CCN_INTERP_PERLIN
} ccnInterpolationMethod;

typedef enum {
	CCN_DIST_MANHATTAN,
	CCN_DIST_EUCLIDEAN,
	CCN_DIST_CHEBYCHEV
} ccnDistanceMethod;

typedef enum {
	CCN_TILE_NOT,
	CCN_TILE_CARTESIAN
} ccnTileMethod;

typedef enum {
	CCN_STORE_SET,
	CCN_STORE_ADD,
	CCN_STORE_SUBTRACT,
	CCN_STORE_MULTIPLY,
	CCN_STORE_DIVIDE
} ccnStoreMethod;

typedef struct {
	uint32_t width, height;
	float *values;
} ccnNoise;

typedef struct {
	ccnTileMethod tileMethod;
	uint32_t xPeriod, yPeriod;
} ccnTileConfiguration;

typedef struct {
	float low, high;
} ccnRange;

typedef struct {
	uint32_t seed;
	int32_t x, y, z;
	ccnTileConfiguration tileConfiguration;
	ccnStoreMethod storeMethod;
	ccnRange range;
} ccnNoiseConfiguration;

// White noise

void ccnGenerateWhiteNoise2D(
	ccnNoise *noise,
	const ccnNoiseConfiguration *configuration);

// Value noise
//
// scale:               the size of a single interpolation interval
// interpolationMethod: the method by which the distance value is interpolated

void ccnGenerateValueNoise1D(
	ccnNoise *noise,
	const ccnNoiseConfiguration *configuration,
	const uint32_t scale,
	const  ccnInterpolationMethod interpolationMethod);

void ccnGenerateValueNoise2D(
	ccnNoise *noise,
	const ccnNoiseConfiguration *configuration,
	const uint32_t scale,
	const ccnInterpolationMethod interpolationMethod);

// Worley noise
//
// points:              the number of points per noise
// n:                   worley noise interpolates to the n-th closest point
// low:                 lowest distance bound for interpolation
// high:                highest distance bound for interpolation
// distanceMethod:      the method by which the distance to a point is calculated
// interpolationMethod: the method by which the distance value is interpolated

void ccnGenerateWorleyNoise2D(
	ccnNoise *noise,
	const ccnNoiseConfiguration *configuration,
	const uint32_t points,
	const uint32_t n,
	const uint32_t low, uint32_t high,
	const ccnDistanceMethod distanceMethod,
	const ccnInterpolationMethod interpolationMethod);

// Perlin noise
//
// scale:               the size of a single interpolation interval
// interpolationMethod: the method by which the distance value is interpolated

void ccnGeneratePerlinNoise1D(
	ccnNoise *noise,
	const ccnNoiseConfiguration *configuration,
	const uint32_t scale,
	const ccnInterpolationMethod interpolationMethod);

void ccnGeneratePerlinNoise2D(
	ccnNoise *noise,
	const ccnNoiseConfiguration *configuration,
	const uint32_t scale,
	const ccnInterpolationMethod interpolationMethod);

#ifdef __cplusplus
}
#endif