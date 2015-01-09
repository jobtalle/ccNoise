//__________________________________________________________________________________//
//                                 _   _       _                                    //
//                                | \ | |     (_)                                   //
//                         ___ ___|  \| | ___  _ ___  ___                           //
//                        / __/ __| . ` |/ _ \| / __|/ _ \                          //
//                       | (_| (__| |\  | (_) | \__ \  __/                          //
//                        \___\___|_| \_|\___/|_|___/\___| 1.0                      //
//                                                                                  //
//              Copyright (C) 2014 - 2015 \ Job Talle (job@ccore.org)               //
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

#define CCN_ERROR_NONE                   0x00
#define CCN_ERROR_INVALID_METHOD         0x01

#define ccnNoiseAllocate(noise, w, h) noise.values = malloc(sizeof(float)* w * h); noise.width = w; noise.height = h;
#define ccnNoiseGet(noise, x, y) noise.values[x + noise.line * y]
#define ccnNoiseSet(noise, x, y, value) noise.values[x + noise.line * y] = value
#define ccnNoiseFree(noise) free(noise.values)

typedef enum {
	CCN_INTERP_LINEAR,
	CCN_INTERP_QUADRATIC,
	CCN_INTERP_QUADRATIC_INVERSE,
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
	unsigned int width, height;
	float *values;
} ccnNoise;

typedef struct {
	ccnTileMethod tileMethod;
	unsigned int xPeriod, yPeriod;
} ccnTileConfiguration;

typedef struct {
	float low, high;
} ccnRange;

typedef struct {
	unsigned int seed;
	int x, y;
	ccnTileConfiguration tileConfiguration;
	ccnStoreMethod storeMethod;
	ccnRange range;
} ccnNoiseConfiguration;

// Create white noise
int ccnGenerateWhiteNoise(
	ccnNoise *noise,
	ccnNoiseConfiguration *configuration);

// Create value noise
int ccnGenerateValueNoise(
	ccnNoise *noise,
	ccnNoiseConfiguration *configuration,
	unsigned int scale,                          // The size of a single interpolation interval
	ccnInterpolationMethod interpolationMethod); // The method by which the distance value is interpolated

// Create worley noise
int ccnGenerateWorleyNoise(
	ccnNoise *noise,
	ccnNoiseConfiguration *configuration,
	unsigned int points,                         // The number of points per noise
	unsigned int n,                              // Worley noise interpolates to the n-th closest point
	int low, int high,                           // Interpolation occurs between the lowest and highest distance
	ccnDistanceMethod distanceMethod,            // The method by which the distance to a point is calculated
	ccnInterpolationMethod interpolationMethod); // The method by which the distance value is interpolated

// Create perlin noise
int ccnGeneratePerlinNoise(
	ccnNoise *noise,
	ccnNoiseConfiguration *configuration,
	unsigned int scale,                          // The size of a single interpolation interval
	ccnInterpolationMethod interpolationMethod); // The method by which the distance value is interpolated

#ifdef __cplusplus
}
#endif