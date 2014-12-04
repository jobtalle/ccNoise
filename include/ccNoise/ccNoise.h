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

#ifdef __cplusplus
extern "C"
{
#endif

// Interpolation methods for noise generation
typedef enum {
	CCN_INTERP_LINEAR,
	CCN_INTERP_COSINE,
	CCN_INTERP_CUBIC
} ccnInterpolationMethod;

typedef enum {
	CCN_DIST_MANHATTAN,
	CCN_DIST_EUCLIDEAN
} ccnDistanceMethod;

// Create an unique as possible ID for a coordinate
unsigned int ccnCoordinateUid(int x, int y);

// Create worley noise
void ccnGenerateWorleyNoise(
	float **buffer,
	unsigned int seed,
	int x, int y,
	unsigned int width, unsigned int height,
	unsigned int points,
	unsigned int n,
	int low, int high,
	float lowValue, float highValue,
	ccnDistanceMethod distanceMethod);

#ifdef __cplusplus
}
#endif