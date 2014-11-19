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

ccRandomizer randomizer;

typedef struct {
	unsigned char r, g, b, a;
} pixelRGBA;

static void generatePerlinNoise()
{
	pixelRGBA *pixels = malloc(sizeof(pixelRGBA)* (WIDTH * HEIGHT));
	float *noise = NULL;
	unsigned int sizes[2] = {WIDTH, HEIGHT};

	//noise = ccnGeneratePerlinNoise(ccrGenerateUint(&randomizer) ^ 42, WIDTH, HEIGHT, 5, 90, 0.5f);
	ccnGenerateValueNoise(ccrGenerateUint(&randomizer) ^ 42, &noise, 2, sizes, 5, 90, 0.5f);

	for(unsigned int i = 0; i < WIDTH * HEIGHT; i++) {
		pixels[i].r = pixels[i].g = pixels[i].b = (unsigned char)(noise[i] * 255.0f);
		pixels[i].a = 255;
	}

	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, WIDTH, HEIGHT, 0, GL_RGBA, GL_UNSIGNED_BYTE, pixels);

	free(pixels);
	free(noise);
}

static void testGridNumberer(void)
{
#define GRIDRADIUS 2

	struct coordinates{
		int x, y, z;
	};

	struct coordinates c;

	c.z = -1;

	for(c.z = -GRIDRADIUS; c.z <= GRIDRADIUS; c.z++) {
		for(c.y = -GRIDRADIUS; c.y <= GRIDRADIUS; c.y++) {
			printf("\n");
			for(c.x = -GRIDRADIUS; c.x <= GRIDRADIUS; c.x++) {
				printf("%d\t", ccnCoordinateUid(3, (int*)&c));
			}
			printf("\n");
		}
		printf("\n");
	}
}

int main(int argc, char **argv)
{
	bool loop = true;
	
	testGridNumberer();

	ccrSeed(&randomizer, (unsigned int)ccTimeNanoseconds());

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
				case CC_KEY_SPACE:
					printf("Rendering perlin noise...\n");
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