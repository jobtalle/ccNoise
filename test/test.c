#include <stdio.h>
#include <math.h>

#include <ccore/display.h>
#include <ccore/window.h>
#include <ccore/opengl.h>
#include <ccore/time.h>

#include <ccNoise/ccNoise.h>
#include <ccRandom/ccRandom.h>

#ifdef WINDOWS
#include <gl/GL.h>
#else
#include <GL/glew.h>
#endif

#define WIDTH  512
#define HEIGHT (256 + 128)

#ifndef _DEBUG
#pragma comment(linker, "/SUBSYSTEM:windows /ENTRY:mainCRTStartup")
#endif

GLuint textureLeftTop, textureRightTop, textureLeftBottom, textureRightBottom;

ccRandomizer32 randomizer;
unsigned int seed;

typedef struct {
	unsigned char r, g, b, a;
} pixelRGBA;

struct osn_context {
	int16_t *perm;
	int16_t *permGradIndex3D;
};

static void generate(int left, int top)
{
	pixelRGBA *pixels = malloc(sizeof(pixelRGBA)* (WIDTH * HEIGHT));

	ccnNoise noise;
	ccnNoiseConfiguration config;

	// Allocate the noise
	ccnNoiseAllocate2D(noise, WIDTH, HEIGHT);

	// Seed & initial store method
	config.seed = seed;
	config.storeMethod = CCN_STORE_SET;

	// Noise coordinates, used for tiling
	config.x = left?1:2;
	config.y = top?7:8;

	// Tile configuration
	config.tileConfiguration.tileMethod = CCN_TILE_CARTESIAN;
	config.tileConfiguration.xPeriod = 2;
	config.tileConfiguration.yPeriod = 2;

	// Noise range
	config.range.low = 0;
	config.range.high = 2;
	
	ccnGenerateOpenSimplex2D(&noise, &config, 32);
	//ccnGenerateValueNoise2D(&noise, &config, 32, CCN_INTERP_CUBIC);

	// Create texture from noise
	unsigned int i;
	for(i = 0; i < WIDTH * HEIGHT; i++) {
		unsigned char v = (unsigned char)(noise.values[i] * 255);

		pixels[i].r = pixels[i].g = pixels[i].b = v;
		//pixels[i].b = pixels[i].r - 160;
	}

	ccnNoiseFree(noise);

	glBindTexture(GL_TEXTURE_2D, left?top?textureLeftTop:textureLeftBottom:top?textureRightTop:textureRightBottom);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, WIDTH, HEIGHT, 0, GL_RGBA, GL_UNSIGNED_BYTE, pixels);

	free(pixels);
}

static void generateLeftTop()
{
	generate(1, 1);
}

static void generateRightTop()
{
	generate(0, 1);
}

static void generateLeftBottom()
{
	generate(1, 0);
}

static void generateRightBottom()
{
	generate(0, 0);
}

int main(int argc, char **argv)
{
	bool loop = true;

	ccrSeed32(&randomizer, (unsigned int)ccTimeNanoseconds());
	seed = ccrGenerateUint32(&randomizer);

	ccDisplayInitialize();

	ccWindowCreate((ccRect){ 0, 0, WIDTH << 1, HEIGHT << 1}, "ccNoise test", CC_WINDOW_FLAG_NORESIZE);
	ccWindowSetCentered();
	//ccWindowSetFullscreen(CC_FULLSCREEN_CURRENT_DISPLAY);

	ccGLContextBind();

	glEnable(GL_TEXTURE_2D);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);

	glGenTextures(1, &textureLeftTop);
	glGenTextures(1, &textureRightTop);
	glGenTextures(1, &textureLeftBottom);
	glGenTextures(1, &textureRightBottom);

	glBindTexture(GL_TEXTURE_2D, textureLeftTop);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glBindTexture(GL_TEXTURE_2D, textureRightTop);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glBindTexture(GL_TEXTURE_2D, textureLeftBottom);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glBindTexture(GL_TEXTURE_2D, textureRightBottom);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

	while(loop) {
		while(ccWindowEventPoll()) {
			if(ccWindowEventGet().type == CC_EVENT_WINDOW_QUIT) loop = false;
			if(ccWindowEventGet().type == CC_EVENT_KEY_DOWN) {
				switch(ccWindowEventGet().keyCode) {
				case CC_KEY_SPACE:
					seed = ccrGenerateUint32(&randomizer);
					printf("Set seed to %d\n", seed);
					break;
				case CC_KEY_1:
					printf("Rendering noise...\n");
					generateLeftTop();
					printf("done.\n");
					break;
				case CC_KEY_2:
					printf("Rendering noise...\n");
					generateRightTop();
					printf("done.\n");
					break;
				case CC_KEY_3:
					printf("Rendering noise...\n");
					generateLeftBottom();
					printf("done.\n");
					break;
				case CC_KEY_4:
					printf("Rendering noise...\n");
					generateRightBottom();
					printf("done.\n");
					break;
				case CC_KEY_5:
					printf("Rendering all...\n");
					seed = ccrGenerateUint32(&randomizer);
					generateLeftTop();
					generateRightTop();
					generateLeftBottom();
					generateRightBottom();
					printf("done.\n");
					break;
				case CC_KEY_ESCAPE:
					loop = false;
					break;
				}
			}
		}

		glClear(GL_COLOR_BUFFER_BIT);

		glBindTexture(GL_TEXTURE_2D, textureLeftTop);
		glBegin(GL_QUADS);
		glTexCoord2f(0.0f, 0.0f);	glVertex2f(-1.0f, 1.0f);
		glTexCoord2f(0.0f, 1.0f);	glVertex2f(-1.0f, 0.0f);
		glTexCoord2f(1.0f, 1.0f);	glVertex2f(0.0f, 0.0f);
		glTexCoord2f(1.0f, 0.0f);	glVertex2f(0.0f, 1.0f);
		glEnd();

		glBindTexture(GL_TEXTURE_2D, textureRightTop);
		glBegin(GL_QUADS);
		glTexCoord2f(0.0f, 0.0f);	glVertex2f(0.0f, 1.0f);
		glTexCoord2f(0.0f, 1.0f);	glVertex2f(0.0f, 0.0f);
		glTexCoord2f(1.0f, 1.0f);	glVertex2f(1.0f, 0.0f);
		glTexCoord2f(1.0f, 0.0f);	glVertex2f(1.0f, 1.0f);
		glEnd();

		glBindTexture(GL_TEXTURE_2D, textureLeftBottom);
		glBegin(GL_QUADS);
		glTexCoord2f(0.0f, 0.0f);	glVertex2f(-1.0f, 0.0f);
		glTexCoord2f(0.0f, 1.0f);	glVertex2f(-1.0f, -1.0f);
		glTexCoord2f(1.0f, 1.0f);	glVertex2f(0.0f, -1.0f);
		glTexCoord2f(1.0f, 0.0f);	glVertex2f(0.0f, 0.0f);
		glEnd();

		glBindTexture(GL_TEXTURE_2D, textureRightBottom);
		glBegin(GL_QUADS);
		glTexCoord2f(0.0f, 0.0f);	glVertex2f(0.0f, 0.0f);
		glTexCoord2f(0.0f, 1.0f);	glVertex2f(0.0f, -1.0f);
		glTexCoord2f(1.0f, 1.0f);	glVertex2f(1.0f, -1.0f);
		glTexCoord2f(1.0f, 0.0f);	glVertex2f(1.0f, 0.0f);
		glEnd();

		ccGLBuffersSwap();
	}

	ccFree();

	return 0;
}
