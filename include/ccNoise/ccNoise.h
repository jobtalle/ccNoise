//__________________________________________________________________________________//
//                                                                                  //
//                                 _   _       _                                    //
//                                | \ | |     (_)                                   //
//                         ___ ___|  \| | ___  _ ___  ___                           //
//                        / __/ __| . ` |/ _ \| / __|/ _ \                          //
//                       | (_| (__| |\  | (_) | \__ \  __/                          //
//                        \___\___|_| \_|\___/|_|___/\___| 1.0                      //
//                                                                                  //
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

#ifdef __cplusplus
extern "C"
{
#endif

// Create an unique as possible ID for a coordinate
unsigned int ccnCoordinateUid(unsigned int dimensions, int *values);

// Creates n-dimensional perlin noise
// persistence 0 - 0.5
void ccnGenerateValueNoise(unsigned int seed, float **buffer, unsigned int dimensions, unsigned int *sizes, unsigned int octaves, unsigned int maxOctave, float octavePersistence);

#ifdef __cplusplus
}
#endif