#include <stdio.h>

#include <ccore/display.h>
#include <ccore/window.h>
#include <ccore/opengl.h>
#include <ccore/time.h>

#include <ccNoise/ccNoise.h>
#include <ccRandom/ccRandom.h>

#include <gl/GL.h>

#define WIDTH  256
#define HEIGHT 256

GLuint textureLeftTop, textureRightTop, textureLeftBottom, textureRightBottom;

ccRandomizer32 randomizer;
unsigned int seed;

typedef struct {
	unsigned char r, g, b, a;
} pixelRGBA;

static void generate(int left, int top)
{
	pixelRGBA *pixels = malloc(sizeof(pixelRGBA)* (WIDTH * HEIGHT));

	ccnNoise noise;
	ccnNoiseConfiguration config;

	ccnNoiseAllocate(noise, WIDTH, HEIGHT);

	config.seed = seed;
	config.range = (ccnRange){ 0.6f, 0.0f };
	config.storeMethod = CCN_STORE_SET;
	config.x = left?0:1;
	config.y = top?0:1;

	config.tileConfiguration.tileMethod = CCN_TILE_CARTESIAN;
	config.tileConfiguration.xPeriod = 2;
	config.tileConfiguration.yPeriod = 2;

	ccnGenerateWorleyNoise(noise, config, 30, 0, 0, 45, CCN_DIST_EUCLIDEAN, CCN_INTERP_COSINE);

	config.range = (ccnRange){ 0.0f, 3.6f };
	config.storeMethod = CCN_STORE_ADD;

	ccnGeneratePerlinNoise(noise, config, 256, CCN_INTERP_PERLIN);

	for(unsigned int i = 0; i < WIDTH * HEIGHT; i++) {
		pixels[i].r = pixels[i].g = pixels[i].b = noise.values[i] > 1.95f?(unsigned char)(noise.values[i] * 255.0f):0;
		pixels[i].a = 255;
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

	ccWindowCreate((ccRect){ 0, 0, WIDTH << 1, HEIGHT << 1 }, "ccNoise test", CC_WINDOW_FLAG_NORESIZE);
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
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glBindTexture(GL_TEXTURE_2D, textureRightTop);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glBindTexture(GL_TEXTURE_2D, textureLeftBottom);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glBindTexture(GL_TEXTURE_2D, textureRightBottom);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

	while(loop) {
		while(ccWindowEventPoll()) {
			if(ccWindowEventGet().type == CC_EVENT_WINDOW_QUIT) loop = false;
			if(ccWindowEventGet().type == CC_EVENT_KEY_DOWN) {
				switch(ccWindowEventGet().keyCode) {
				case CC_KEY_SPACE:
					seed = ccrGenerateUint32(&randomizer) ^ 42;
					printf("Set seed to %d\n", seed);
					break;
				case CC_KEY_1:
					printf("Rendering worley noise...\n");
					generateLeftTop();
					printf("done.\n");
					break;
				case CC_KEY_2:
					printf("Rendering worley noise...\n");
					generateRightTop();
					printf("done.\n");
					break;
				case CC_KEY_3:
					printf("Rendering worley noise...\n");
					generateLeftBottom();
					printf("done.\n");
					break;
				case CC_KEY_4:
					printf("Rendering worley noise...\n");
					generateRightBottom();
					printf("done.\n");
					break;
				case CC_KEY_5:
					printf("Rendering all...\n");
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