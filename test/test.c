#include <stdio.h>

#include <ccore/window.h>
#include <ccore/display.h>
#include <ccore/opengl.h>

#include <ccNoise/ccNoise.h>

int main(int argc, char **argv)
{
	bool loop = true;

	ccDisplayInitialize();

	ccWindowCreate((ccRect){ 0, 0, 800, 800 }, "ccNoise test", CC_WINDOW_FLAG_NORESIZE);
	ccWindowSetCentered();

	ccGLContextBind();

	while(loop) {
		while(ccWindowEventPoll()) {
			if(ccWindowEventGet().type == CC_EVENT_WINDOW_QUIT) loop = false;
		}

		ccGLBuffersSwap();
	}

	ccFree();
}