#include <ccRandom/ccRandom.h>

#define CCR_MULTIPLY_32           69069
#define CCR_ADD_32                1
#define CCR_MULTIPLY_64           6364136223846793005
#define CCR_ADD_64                1
#define CCR_MULTIPLY_COORDINATE_A 134775813
#define CCR_MULTIPLY_COORDINATE_B 1103515245

#define CCR_LCG32 *randomizer = *randomizer * CCR_MULTIPLY_32 + CCR_ADD_32
#define CCR_LCG64 *randomizer = *randomizer * CCR_MULTIPLY_64 + CCR_ADD_64
#define CCR_CRG (((x ^ y) * CCR_MULTIPLY_COORDINATE_A) ^ (seed + x)) * (((CCR_MULTIPLY_COORDINATE_B * x) << 16) ^ (CCR_MULTIPLY_COORDINATE_B * y) - CCR_MULTIPLY_COORDINATE_A)

void ccrSeed32(ccRandomizer32 *randomizer, uint32_t seed)
{
	*randomizer = seed;
}

void ccrSeed64(ccRandomizer64 *randomizer, uint64_t seed)
{
	*randomizer = seed;
}

uint32_t ccrGenerateUint32(ccRandomizer32 *randomizer)
{
	return CCR_LCG32;
}

uint64_t ccrGenerateUint64(ccRandomizer64 *randomizer)
{
	return CCR_LCG64;
}

uint32_t ccrGenerateUintCoordinate(uint32_t seed, int32_t x, int32_t y)
{
	return CCR_CRG;
}

float ccrGenerateFloat32(ccRandomizer32 *randomizer)
{
	return (float)(CCR_LCG32) / UINT32_MAX;
}

float ccrGenerateFloat64(ccRandomizer64 *randomizer)
{
	return (float)(CCR_LCG64) / UINT64_MAX;
}

float ccrGenerateFloatCoordinate(uint32_t seed, int32_t x, int32_t y)
{
	return (float)(CCR_CRG) / UINT32_MAX;
}

double ccrGenerateDouble32(ccRandomizer32 *randomizer)
{
	return (double)(CCR_LCG32) / UINT32_MAX;
}

double ccrGenerateDouble64(ccRandomizer64 *randomizer)
{
	return (double)(CCR_LCG64) / UINT64_MAX;
}

double ccrGenerateDoubleCoordinate(uint32_t seed, int32_t x, int32_t y)
{
	return (double)(CCR_CRG) / UINT32_MAX;
}