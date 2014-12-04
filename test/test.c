#include <stdio.h>

#include <ccore/display.h>
#include <ccore/window.h>
#include <ccore/opengl.h>
#include <ccore/time.h>

#include <ccNoise/ccNoise.h>
#include <ccRandom/ccRandom.h>

#include <gl/GL.h>

#define WIDTH  512
#define HEIGHT 512

GLuint texture;

ccRandomizer32 randomizer;

typedef struct {
	unsigned char r, g, b, a;
} pixelRGBA;

static void generatePerlinNoise()
{
	pixelRGBA *pixels = malloc(sizeof(pixelRGBA)* (WIDTH * HEIGHT));
	float *noise = NULL;

	ccnGenerateWorleyNoise(&noise, ccrGenerateUint32(&randomizer) ^ 42, 0, 0, WIDTH, HEIGHT, 24, 0, 0, 140, 0.1f, 1.0f, CCN_DIST_EUCLIDEAN);

	for(unsigned int i = 0; i < WIDTH * HEIGHT; i++) {
		pixels[i].r = pixels[i].g = pixels[i].b = (unsigned char)(noise[i] * 255.0f);
		pixels[i].a = 255;
	}

	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, WIDTH, HEIGHT, 0, GL_RGBA, GL_UNSIGNED_BYTE, pixels);

	free(pixels);
	free(noise);
}

int main(int argc, char **argv)
{
	bool loop = true;

	ccrSeed32(&randomizer, (unsigned int)ccTimeNanoseconds());

	ccDisplayInitialize();

	ccWindowCreate((ccRect){ 0, 0, WIDTH, HEIGHT }, "ccNoise test", CC_WINDOW_FLAG_NORESIZE);
	ccWindowSetCentered();

	ccGLContextBind();

	glEnable(GL_TEXTURE_2D);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);

	glGenTextures(1, &texture);
	glBindTexture(GL_TEXTURE_2D, texture);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

	while(loop) {
		while(ccWindowEventPoll()) {
			if(ccWindowEventGet().type == CC_EVENT_WINDOW_QUIT) loop = false;
			if(ccWindowEventGet().type == CC_EVENT_KEY_DOWN) {
				switch(ccWindowEventGet().keyCode) {
				case CC_KEY_1:
					printf("Rendering worley noise...\n");
					generatePerlinNoise();
					printf("done.\n");
					break;
				}
			}
		}

		glClear(GL_COLOR_BUFFER_BIT);

		glBegin(GL_QUADS);
		glTexCoord2f(0.0f, 0.0f);	glVertex2f(-1.0f, 1.0f);
		glTexCoord2f(0.0f, 1.0f);	glVertex2f(-1.0f, -1.0f);
		glTexCoord2f(1.0f, 1.0f);	glVertex2f(1.0f, -1.0f);
		glTexCoord2f(1.0f, 0.0f);	glVertex2f(1.0f, 1.0f);
		glEnd();

		ccGLBuffersSwap();
	}

	ccFree();
}