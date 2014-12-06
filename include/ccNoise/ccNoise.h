//__________________________________________________________________________________//
//                                 _   _       _                                    //
//                                | \ | |     (_)                                   //
//                         ___ ___|  \| | ___  _ ___  ___                           //
//                        / __/ __| . ` |/ _ \| / __|/ _ \                          //
//                       | (_| (__| |\  | (_) | \__ \  __/                          //
//                        \___\___|_| \_|\___/|_|___/\___| 1.0                      //
//                                                                                  //
//             Copyright (C) 2014 \ Job Talle (job@ccore.org)                       //
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

#include <ccore/types.h>

#ifdef __cplusplus
extern "C"
{
#endif

typedef enum {
	CCN_INTERP_LINEAR,
	CCN_INTERP_QUADRATIC,
	CCN_INTERP_QUADRATIC_INVERSE,
	CCN_INTERP_COSINE
} ccnInterpolationMethod;

typedef enum {
	CCN_DIST_MANHATTAN,
	CCN_DIST_EUCLIDEAN,
	CCN_DIST_CHEBYCHEV
} ccnDistanceMethod;

// Create an unique as possible ID for a coordinate
unsigned int ccnCoordinateUid(int x, int y);

// Create worley noise
void ccnGenerateWorleyNoise(
	float **buffer,                              // The buffer to store the generated values in
	unsigned int seed,                           // The random seed
	int x, int y,                                // Adjecent coordinates will tile seamlessly
	unsigned int width, unsigned int height,     // Noise dimensions
	unsigned int points,                         // The number of points per noise
	unsigned int n,                              // Worley noise interpolates to the n-th closest point
	int low, int high,                           // Interpolation occurs between the lowest and highest distance
	float lowValue, float highValue,             // Minimum and maximum value of the noise
	ccnDistanceMethod distanceMethod,            // The method by which the distance to a point is calculated
	ccnInterpolationMethod interpolationMethod); // The method by which the distance value is interpolated

// Create value noise
void ccnGenerateFractalNoise(
	float **buffer,                              // The buffer to store the generated values in
	unsigned int seed,                           // The random seed
	bool makeTileable,                           // Make noise tileable, costs some overhead
	int x, int y,                                // Adjecent coordinates will tile seamlessly
	unsigned int width, unsigned int height,     // Noise dimensions
	unsigned int octaves,                        // The number of times to add noises to the noise, 2 log maxOctaves for maximum detail
	unsigned int maxOctave,                      // The largest interpolation distance, halved for each octave
	float persistence,                           // The factor below to multiply each subsequent octave's influence with
	ccnInterpolationMethod interpolationMethod); // The method by which the distance value is interpolated

#ifdef __cplusplus
}
#endif