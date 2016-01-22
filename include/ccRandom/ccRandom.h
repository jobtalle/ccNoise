//__________________________________________________________________________________//
//                                                                                  //
//                       _____                 _                                    //
//                      |  __ \               | |                                   //
//               ___ ___| |__) |__ _ _ __   __| | ___  _ __ ___                     //
//              / __/ __|  _  // _` | '_ \ / _` |/ _ \| '_ ` _ \                    //
//             | (_| (__| | \ \ (_| | | | | (_| | (_) | | | | | |                   //
//              \___\___|_|  \_\__,_|_| |_|\__,_|\___/|_| |_| |_| 1.0               //
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

typedef uint32_t ccRandomizer32;
typedef uint64_t ccRandomizer64;

void ccrSeed32(ccRandomizer32 *randomizer, uint32_t seed);
void ccrSeed64(ccRandomizer64 *randomizer, uint64_t seed);

uint32_t ccrGenerateUint32(ccRandomizer32 *randomizer);
uint64_t ccrGenerateUint64(ccRandomizer64 *randomizer);
uint32_t ccrGenerateUintCoordinate(uint32_t seed, int32_t x, int32_t y);
float ccrGenerateFloat32(ccRandomizer32 *randomizer);
float ccrGenerateFloat64(ccRandomizer64 *randomizer);
float ccrGenerateFloatCoordinate(uint32_t seed, int32_t x, int32_t y);
double ccrGenerateDouble32(ccRandomizer32 *randomizer);
double ccrGenerateDouble64(ccRandomizer64 *randomizer);
double ccrGenerateDoubleCoordinate(uint32_t seed, int32_t x, int32_t y);

#ifdef __cplusplus
}
#endif