#include <ccTrigonometry/ccTrigonometry.h>

// Clamp

void ccTriClampRad(float *value)
{
	for(; *value < -CC_TRI_PI_F; *value += CC_TRI_PI_DOUBLE_F);
	for(; *value > CC_TRI_PI_F; *value -= CC_TRI_PI_DOUBLE_F);
}

void ccTriClampf(float *value, float low, float high)
{
	if(*value < low) {
		*value = low;
	}
	else if(*value > high) {
		*value = high;
	}
}

// Delta

float ccTriDeltaRad(float a, float b)
{
	float delta = a - b;

	ccTriClampRad(&delta);

	return delta;
}

// Lerp

// Perp

void ccTriPerpRad(float *value, float target, float factor)
{
	*value += ccTriDeltaRad(target, *value) * factor;
}

void ccTriPerpf(float *value, float target, float factor)
{
	*value += (target - *value) * factor;
}