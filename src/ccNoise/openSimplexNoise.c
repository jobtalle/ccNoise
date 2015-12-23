#include "ccNoiseInternal.h"

/*
* OpenSimplex (Simplectic) Noise in C.
* Ported by Stephen M. Cameron from Kurt Spencer's java implementation
*
* v1.1 (October 5, 2014)
* - Added 2D and 4D implementations.
* - Proper gradient sets for all dimensions, from a
*   dimensionally-generalizable scheme with an actual
*   rhyme and reason behind it.
* - Removed default permutation array in favor of
*   default seed.
* - Changed seed-based constructor to be independent
*   of any particular randomization library, so results
*   will be the same when ported to other languages.
*/
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>

#ifndef inline
#define INLINE _inline
#else
#define INLINE inline
#endif

struct osn_context {
	int16_t *perm;
	int16_t *permGradIndex3D;
};

int open_simplex_noise(int64_t seed, struct osn_context **ctx);
void open_simplex_noise_free(struct osn_context *ctx);
float open_simplex_noise2(struct osn_context *ctx, float x, float y);
float open_simplex_noise3(struct osn_context *ctx, float x, float y, float z);
float open_simplex_noise4(struct osn_context *ctx, float x, float y, float z, float w);

#define STRETCH_CONSTANT_2D -0.211324865405187     /* (1 / sqrt(2 + 1) - 1 ) / 2; */
#define SQUISH_CONSTANT_2D   0.366025403784439     /* (sqrt(2 + 1) -1) / 2; */
#define STRETCH_CONSTANT_3D -1.0 / 6.0             /* (1 / sqrt(3 + 1) - 1) / 3; */
#define SQUISH_CONSTANT_3D   1.0 / 3.0             /* (sqrt(3+1)-1)/3; */

#define NORM_CONSTANT_2D 81.406388f                /* 94 * (sqrt(3) / 2) */
#define NORM_CONSTANT_3D 103.0

#define DEFAULT_SEED (0LL)

#define ARRAYSIZE(x) (sizeof((x)) / sizeof((x)[0]))

/* 
 * Gradients for 2D. They approximate the directions to the
 * vertices of an octagon from the center.
 */
static const int8_t gradients2D[] = {
	 5,  2,    2,  5,
	-5,  2,   -2,  5,
	 5, -2,    2, -5,
	-5, -2,   -2, -5,
};

/*	
 * Gradients for 3D. They approximate the directions to the
 * vertices of a rhombicuboctahedron from the center, skewed so
 * that the triangular and square facets can be inscribed inside
 * circles of the same radius.
 */
static const signed char gradients3D[] = {
	-11,  4,  4,     -4,  11,  4,    -4,  4,  11,
	 11,  4,  4,      4,  11,  4,     4,  4,  11,
	-11, -4,  4,     -4, -11,  4,    -4, -4,  11,
	 11, -4,  4,      4, -11,  4,     4, -4,  11,
	-11,  4, -4,     -4,  11, -4,    -4,  4, -11,
	 11,  4, -4,      4,  11, -4,     4,  4, -11,
	-11, -4, -4,     -4, -11, -4,    -4, -4, -11,
	 11, -4, -4,      4, -11, -4,     4, -4, -11,
};

static float extrapolate2(struct osn_context *ctx, int xsb, int ysb, float dx, float dy)
{
	int16_t *perm = ctx->perm;
	int index = perm[(perm[xsb & 0xFF] + ysb) & 0xFF] & 0x0E;
	return gradients2D[index] * dx
		+ gradients2D[index + 1] * dy;
}

static float extrapolate3(struct osn_context *ctx, int xsb, int ysb, int zsb, float dx, float dy, float dz)
{
	int16_t *perm = ctx->perm;
	int16_t *permGradIndex3D = ctx->permGradIndex3D;
	int index = permGradIndex3D[(perm[(perm[xsb & 0xFF] + ysb) & 0xFF] + zsb) & 0xFF];
	return gradients3D[index] * dx
		+ gradients3D[index + 1] * dy
		+ gradients3D[index + 2] * dz;
}

static INLINE int fastFloor(float x) {
	int xi = (int)x;
	return x < xi?xi - 1:xi;
}

static int allocate_perm(struct osn_context *ctx, int nperm, int ngrad)
{
	if(ctx->perm)
		free(ctx->perm);
	if(ctx->permGradIndex3D)
		free(ctx->permGradIndex3D);
	ctx->perm = (int16_t *)malloc(sizeof(*ctx->perm) * nperm);
	if(!ctx->perm)
		return -ENOMEM;
	ctx->permGradIndex3D = (int16_t *)malloc(sizeof(*ctx->permGradIndex3D) * ngrad);
	if(!ctx->permGradIndex3D) {
		free(ctx->perm);
		return -ENOMEM;
	}
	return 0;
}

int open_simplex_noise_init_perm(struct osn_context *ctx, int16_t p[], int nelements)
{
	int i, rc;

	rc = allocate_perm(ctx, nelements, 256);
	if(rc)
		return rc;
	memcpy(ctx->perm, p, sizeof(*ctx->perm) * nelements);

	for(i = 0; i < 256; i++) {
		/* Since 3D has 24 gradients, simple bitmask won't work, so precompute modulo array. */
		ctx->permGradIndex3D[i] = (int16_t)((ctx->perm[i] % (ARRAYSIZE(gradients3D) / 3)) * 3);
	}
	return 0;
}

/*
* Initializes using a permutation array generated from a 64-bit seed.
* Generates a proper permutation (i.e. doesn't merely perform N successive pair
* swaps on a base array).  Uses a simple 64-bit LCG.
*/
int open_simplex_noise(int64_t seed, struct osn_context **ctx)
{
	int rc;
	int16_t source[256];
	int i;
	int16_t *perm;
	int16_t *permGradIndex3D;
	int r;

	*ctx = (struct osn_context *) malloc(sizeof(**ctx));
	if(!(*ctx))
		return -ENOMEM;
	(*ctx)->perm = NULL;
	(*ctx)->permGradIndex3D = NULL;

	rc = allocate_perm(*ctx, 256, 256);
	if(rc) {
		free(*ctx);
		return rc;
	}

	perm = (*ctx)->perm;
	permGradIndex3D = (*ctx)->permGradIndex3D;

	for(i = 0; i < 256; i++)
		source[i] = (int16_t)i;
	seed = seed * 6364136223846793005LL + 1442695040888963407LL;
	seed = seed * 6364136223846793005LL + 1442695040888963407LL;
	seed = seed * 6364136223846793005LL + 1442695040888963407LL;
	for(i = 255; i >= 0; i--) {
		seed = seed * 6364136223846793005LL + 1442695040888963407LL;
		r = (int)((seed + 31) % (i + 1));
		if(r < 0)
			r += (i + 1);
		perm[i] = source[r];
		permGradIndex3D[i] = (short)((perm[i] % (ARRAYSIZE(gradients3D) / 3)) * 3);
		source[r] = source[i];
	}
	return 0;
}

void open_simplex_noise_free(struct osn_context *ctx)
{
	if(!ctx)
		return;
	if(ctx->perm) {
		free(ctx->perm);
		ctx->perm = NULL;
	}
	if(ctx->permGradIndex3D) {
		free(ctx->permGradIndex3D);
		ctx->permGradIndex3D = NULL;
	}
	free(ctx);
}

/* 2D OpenSimplex (Simplectic) Noise. */
float open_simplex_noise2(struct osn_context *ctx, float x, float y)
{

	/* Place input coordinates onto grid. */
	float stretchOffset = (x + y) * STRETCH_CONSTANT_2D;
	float xs = x + stretchOffset;
	float ys = y + stretchOffset;

	/* Floor to get grid coordinates of rhombus (stretched square) super-cell origin. */
	int xsb = fastFloor(xs);
	int ysb = fastFloor(ys);

	/* Skew out to get actual coordinates of rhombus origin. We'll need these later. */
	float squishOffset = (xsb + ysb) * SQUISH_CONSTANT_2D;
	float xb = xsb + squishOffset;
	float yb = ysb + squishOffset;

	/* Compute grid coordinates relative to rhombus origin. */
	float xins = xs - xsb;
	float yins = ys - ysb;

	/* Sum those together to get a value that determines which region we're in. */
	float inSum = xins + yins;

	/* Positions relative to origin point. */
	float dx0 = x - xb;
	float dy0 = y - yb;

	/* We'll be defining these inside the next block and using them afterwards. */
	float dx_ext, dy_ext;
	int xsv_ext, ysv_ext;

	float dx1;
	float dy1;
	float attn1;
	float dx2;
	float dy2;
	float attn2;
	float zins;
	float attn0;
	float attn_ext;

	float value = 0;

	/* Contribution (1,0) */
	dx1 = dx0 - 1 - SQUISH_CONSTANT_2D;
	dy1 = dy0 - 0 - SQUISH_CONSTANT_2D;
	attn1 = 2 - dx1 * dx1 - dy1 * dy1;
	if(attn1 > 0) {
		attn1 *= attn1;
		value += attn1 * attn1 * extrapolate2(ctx, xsb + 1, ysb + 0, dx1, dy1);
	}

	/* Contribution (0,1) */
	dx2 = dx0 - 0 - SQUISH_CONSTANT_2D;
	dy2 = dy0 - 1 - SQUISH_CONSTANT_2D;
	attn2 = 2 - dx2 * dx2 - dy2 * dy2;
	if(attn2 > 0) {
		attn2 *= attn2;
		value += attn2 * attn2 * extrapolate2(ctx, xsb + 0, ysb + 1, dx2, dy2);
	}

	if(inSum <= 1) { /* We're inside the triangle (2-Simplex) at (0,0) */
		zins = 1 - inSum;
		if(zins > xins || zins > yins) { /* (0,0) is one of the closest two triangular vertices */
			if(xins > yins) {
				xsv_ext = xsb + 1;
				ysv_ext = ysb - 1;
				dx_ext = dx0 - 1;
				dy_ext = dy0 + 1;
			}
			else {
				xsv_ext = xsb - 1;
				ysv_ext = ysb + 1;
				dx_ext = dx0 + 1;
				dy_ext = dy0 - 1;
			}
		}
		else { /* (1,0) and (0,1) are the closest two vertices. */
			xsv_ext = xsb + 1;
			ysv_ext = ysb + 1;
			dx_ext = dx0 - 1 - 2 * SQUISH_CONSTANT_2D;
			dy_ext = dy0 - 1 - 2 * SQUISH_CONSTANT_2D;
		}
	}
	else { /* We're inside the triangle (2-Simplex) at (1,1) */
		zins = 2 - inSum;
		if(zins < xins || zins < yins) { /* (0,0) is one of the closest two triangular vertices */
			if(xins > yins) {
				xsv_ext = xsb + 2;
				ysv_ext = ysb + 0;
				dx_ext = dx0 - 2 - 2 * SQUISH_CONSTANT_2D;
				dy_ext = dy0 + 0 - 2 * SQUISH_CONSTANT_2D;
			}
			else {
				xsv_ext = xsb + 0;
				ysv_ext = ysb + 2;
				dx_ext = dx0 + 0 - 2 * SQUISH_CONSTANT_2D;
				dy_ext = dy0 - 2 - 2 * SQUISH_CONSTANT_2D;
			}
		}
		else { /* (1,0) and (0,1) are the closest two vertices. */
			dx_ext = dx0;
			dy_ext = dy0;
			xsv_ext = xsb;
			ysv_ext = ysb;
		}
		xsb += 1;
		ysb += 1;
		dx0 = dx0 - 1 - 2 * SQUISH_CONSTANT_2D;
		dy0 = dy0 - 1 - 2 * SQUISH_CONSTANT_2D;
	}

	/* Contribution (0,0) or (1,1) */
	attn0 = 2 - dx0 * dx0 - dy0 * dy0;
	if(attn0 > 0) {
		attn0 *= attn0;
		value += attn0 * attn0 * extrapolate2(ctx, xsb, ysb, dx0, dy0);
	}

	/* Extra Vertex */
	attn_ext = 2 - dx_ext * dx_ext - dy_ext * dy_ext;
	if(attn_ext > 0) {
		attn_ext *= attn_ext;
		value += attn_ext * attn_ext * extrapolate2(ctx, xsv_ext, ysv_ext, dx_ext, dy_ext);
	}

	return value / NORM_CONSTANT_2D + 0.5f;
}

/*
* 3D OpenSimplex (Simplectic) Noise
*/
float open_simplex_noise3(struct osn_context *ctx, float x, float y, float z)
{

	/* Place input coordinates on simplectic honeycomb. */
	float stretchOffset = (x + y + z) * STRETCH_CONSTANT_3D;
	float xs = x + stretchOffset;
	float ys = y + stretchOffset;
	float zs = z + stretchOffset;

	/* Floor to get simplectic honeycomb coordinates of rhombohedron (stretched cube) super-cell origin. */
	int xsb = fastFloor(xs);
	int ysb = fastFloor(ys);
	int zsb = fastFloor(zs);

	/* Skew out to get actual coordinates of rhombohedron origin. We'll need these later. */
	float squishOffset = (xsb + ysb + zsb) * SQUISH_CONSTANT_3D;
	float xb = xsb + squishOffset;
	float yb = ysb + squishOffset;
	float zb = zsb + squishOffset;

	/* Compute simplectic honeycomb coordinates relative to rhombohedral origin. */
	float xins = xs - xsb;
	float yins = ys - ysb;
	float zins = zs - zsb;

	/* Sum those together to get a value that determines which region we're in. */
	float inSum = xins + yins + zins;

	/* Positions relative to origin point. */
	float dx0 = x - xb;
	float dy0 = y - yb;
	float dz0 = z - zb;

	/* We'll be defining these inside the next block and using them afterwards. */
	float dx_ext0, dy_ext0, dz_ext0;
	float dx_ext1, dy_ext1, dz_ext1;
	int xsv_ext0, ysv_ext0, zsv_ext0;
	int xsv_ext1, ysv_ext1, zsv_ext1;

	float wins;
	int8_t c, c1, c2;
	int8_t aPoint, bPoint;
	float aScore, bScore;
	int aIsFurtherSide;
	int bIsFurtherSide;
	float p1, p2, p3;
	float score;
	float attn0, attn1, attn2, attn3, attn4, attn5, attn6;
	float dx1, dy1, dz1;
	float dx2, dy2, dz2;
	float dx3, dy3, dz3;
	float dx4, dy4, dz4;
	float dx5, dy5, dz5;
	float dx6, dy6, dz6;
	float attn_ext0, attn_ext1;

	float value = 0;
	if(inSum <= 1) { /* We're inside the tetrahedron (3-Simplex) at (0,0,0) */

		/* Determine which two of (0,0,1), (0,1,0), (1,0,0) are closest. */
		aPoint = 0x01;
		aScore = xins;
		bPoint = 0x02;
		bScore = yins;
		if(aScore >= bScore && zins > bScore) {
			bScore = zins;
			bPoint = 0x04;
		}
		else if(aScore < bScore && zins > aScore) {
			aScore = zins;
			aPoint = 0x04;
		}

		/* Now we determine the two lattice points not part of the tetrahedron that may contribute.
		This depends on the closest two tetrahedral vertices, including (0,0,0) */
		wins = 1 - inSum;
		if(wins > aScore || wins > bScore) { /* (0,0,0) is one of the closest two tetrahedral vertices. */
			c = (bScore > aScore?bPoint:aPoint); /* Our other closest vertex is the closest out of a and b. */

			if((c & 0x01) == 0) {
				xsv_ext0 = xsb - 1;
				xsv_ext1 = xsb;
				dx_ext0 = dx0 + 1;
				dx_ext1 = dx0;
			}
			else {
				xsv_ext0 = xsv_ext1 = xsb + 1;
				dx_ext0 = dx_ext1 = dx0 - 1;
			}

			if((c & 0x02) == 0) {
				ysv_ext0 = ysv_ext1 = ysb;
				dy_ext0 = dy_ext1 = dy0;
				if((c & 0x01) == 0) {
					ysv_ext1 -= 1;
					dy_ext1 += 1;
				}
				else {
					ysv_ext0 -= 1;
					dy_ext0 += 1;
				}
			}
			else {
				ysv_ext0 = ysv_ext1 = ysb + 1;
				dy_ext0 = dy_ext1 = dy0 - 1;
			}

			if((c & 0x04) == 0) {
				zsv_ext0 = zsb;
				zsv_ext1 = zsb - 1;
				dz_ext0 = dz0;
				dz_ext1 = dz0 + 1;
			}
			else {
				zsv_ext0 = zsv_ext1 = zsb + 1;
				dz_ext0 = dz_ext1 = dz0 - 1;
			}
		}
		else { /* (0,0,0) is not one of the closest two tetrahedral vertices. */
			c = (int8_t)(aPoint | bPoint); /* Our two extra vertices are determined by the closest two. */

			if((c & 0x01) == 0) {
				xsv_ext0 = xsb;
				xsv_ext1 = xsb - 1;
				dx_ext0 = dx0 - 2 * SQUISH_CONSTANT_3D;
				dx_ext1 = dx0 + 1 - SQUISH_CONSTANT_3D;
			}
			else {
				xsv_ext0 = xsv_ext1 = xsb + 1;
				dx_ext0 = dx0 - 1 - 2 * SQUISH_CONSTANT_3D;
				dx_ext1 = dx0 - 1 - SQUISH_CONSTANT_3D;
			}

			if((c & 0x02) == 0) {
				ysv_ext0 = ysb;
				ysv_ext1 = ysb - 1;
				dy_ext0 = dy0 - 2 * SQUISH_CONSTANT_3D;
				dy_ext1 = dy0 + 1 - SQUISH_CONSTANT_3D;
			}
			else {
				ysv_ext0 = ysv_ext1 = ysb + 1;
				dy_ext0 = dy0 - 1 - 2 * SQUISH_CONSTANT_3D;
				dy_ext1 = dy0 - 1 - SQUISH_CONSTANT_3D;
			}

			if((c & 0x04) == 0) {
				zsv_ext0 = zsb;
				zsv_ext1 = zsb - 1;
				dz_ext0 = dz0 - 2 * SQUISH_CONSTANT_3D;
				dz_ext1 = dz0 + 1 - SQUISH_CONSTANT_3D;
			}
			else {
				zsv_ext0 = zsv_ext1 = zsb + 1;
				dz_ext0 = dz0 - 1 - 2 * SQUISH_CONSTANT_3D;
				dz_ext1 = dz0 - 1 - SQUISH_CONSTANT_3D;
			}
		}

		/* Contribution (0,0,0) */
		attn0 = 2 - dx0 * dx0 - dy0 * dy0 - dz0 * dz0;
		if(attn0 > 0) {
			attn0 *= attn0;
			value += attn0 * attn0 * extrapolate3(ctx, xsb + 0, ysb + 0, zsb + 0, dx0, dy0, dz0);
		}

		/* Contribution (1,0,0) */
		dx1 = dx0 - 1 - SQUISH_CONSTANT_3D;
		dy1 = dy0 - 0 - SQUISH_CONSTANT_3D;
		dz1 = dz0 - 0 - SQUISH_CONSTANT_3D;
		attn1 = 2 - dx1 * dx1 - dy1 * dy1 - dz1 * dz1;
		if(attn1 > 0) {
			attn1 *= attn1;
			value += attn1 * attn1 * extrapolate3(ctx, xsb + 1, ysb + 0, zsb + 0, dx1, dy1, dz1);
		}

		/* Contribution (0,1,0) */
		dx2 = dx0 - 0 - SQUISH_CONSTANT_3D;
		dy2 = dy0 - 1 - SQUISH_CONSTANT_3D;
		dz2 = dz1;
		attn2 = 2 - dx2 * dx2 - dy2 * dy2 - dz2 * dz2;
		if(attn2 > 0) {
			attn2 *= attn2;
			value += attn2 * attn2 * extrapolate3(ctx, xsb + 0, ysb + 1, zsb + 0, dx2, dy2, dz2);
		}

		/* Contribution (0,0,1) */
		dx3 = dx2;
		dy3 = dy1;
		dz3 = dz0 - 1 - SQUISH_CONSTANT_3D;
		attn3 = 2 - dx3 * dx3 - dy3 * dy3 - dz3 * dz3;
		if(attn3 > 0) {
			attn3 *= attn3;
			value += attn3 * attn3 * extrapolate3(ctx, xsb + 0, ysb + 0, zsb + 1, dx3, dy3, dz3);
		}
	}
	else if(inSum >= 2) { /* We're inside the tetrahedron (3-Simplex) at (1,1,1) */

		/* Determine which two tetrahedral vertices are the closest, out of (1,1,0), (1,0,1), (0,1,1) but not (1,1,1). */
		aPoint = 0x06;
		aScore = xins;
		bPoint = 0x05;
		bScore = yins;
		if(aScore <= bScore && zins < bScore) {
			bScore = zins;
			bPoint = 0x03;
		}
		else if(aScore > bScore && zins < aScore) {
			aScore = zins;
			aPoint = 0x03;
		}

		/* Now we determine the two lattice points not part of the tetrahedron that may contribute.
		This depends on the closest two tetrahedral vertices, including (1,1,1) */
		wins = 3 - inSum;
		if(wins < aScore || wins < bScore) { /* (1,1,1) is one of the closest two tetrahedral vertices. */
			c = (bScore < aScore?bPoint:aPoint); /* Our other closest vertex is the closest out of a and b. */

			if((c & 0x01) != 0) {
				xsv_ext0 = xsb + 2;
				xsv_ext1 = xsb + 1;
				dx_ext0 = dx0 - 2 - 3 * SQUISH_CONSTANT_3D;
				dx_ext1 = dx0 - 1 - 3 * SQUISH_CONSTANT_3D;
			}
			else {
				xsv_ext0 = xsv_ext1 = xsb;
				dx_ext0 = dx_ext1 = dx0 - 3 * SQUISH_CONSTANT_3D;
			}

			if((c & 0x02) != 0) {
				ysv_ext0 = ysv_ext1 = ysb + 1;
				dy_ext0 = dy_ext1 = dy0 - 1 - 3 * SQUISH_CONSTANT_3D;
				if((c & 0x01) != 0) {
					ysv_ext1 += 1;
					dy_ext1 -= 1;
				}
				else {
					ysv_ext0 += 1;
					dy_ext0 -= 1;
				}
			}
			else {
				ysv_ext0 = ysv_ext1 = ysb;
				dy_ext0 = dy_ext1 = dy0 - 3 * SQUISH_CONSTANT_3D;
			}

			if((c & 0x04) != 0) {
				zsv_ext0 = zsb + 1;
				zsv_ext1 = zsb + 2;
				dz_ext0 = dz0 - 1 - 3 * SQUISH_CONSTANT_3D;
				dz_ext1 = dz0 - 2 - 3 * SQUISH_CONSTANT_3D;
			}
			else {
				zsv_ext0 = zsv_ext1 = zsb;
				dz_ext0 = dz_ext1 = dz0 - 3 * SQUISH_CONSTANT_3D;
			}
		}
		else { /* (1,1,1) is not one of the closest two tetrahedral vertices. */
			c = (int8_t)(aPoint & bPoint); /* Our two extra vertices are determined by the closest two. */

			if((c & 0x01) != 0) {
				xsv_ext0 = xsb + 1;
				xsv_ext1 = xsb + 2;
				dx_ext0 = dx0 - 1 - SQUISH_CONSTANT_3D;
				dx_ext1 = dx0 - 2 - 2 * SQUISH_CONSTANT_3D;
			}
			else {
				xsv_ext0 = xsv_ext1 = xsb;
				dx_ext0 = dx0 - SQUISH_CONSTANT_3D;
				dx_ext1 = dx0 - 2 * SQUISH_CONSTANT_3D;
			}

			if((c & 0x02) != 0) {
				ysv_ext0 = ysb + 1;
				ysv_ext1 = ysb + 2;
				dy_ext0 = dy0 - 1 - SQUISH_CONSTANT_3D;
				dy_ext1 = dy0 - 2 - 2 * SQUISH_CONSTANT_3D;
			}
			else {
				ysv_ext0 = ysv_ext1 = ysb;
				dy_ext0 = dy0 - SQUISH_CONSTANT_3D;
				dy_ext1 = dy0 - 2 * SQUISH_CONSTANT_3D;
			}

			if((c & 0x04) != 0) {
				zsv_ext0 = zsb + 1;
				zsv_ext1 = zsb + 2;
				dz_ext0 = dz0 - 1 - SQUISH_CONSTANT_3D;
				dz_ext1 = dz0 - 2 - 2 * SQUISH_CONSTANT_3D;
			}
			else {
				zsv_ext0 = zsv_ext1 = zsb;
				dz_ext0 = dz0 - SQUISH_CONSTANT_3D;
				dz_ext1 = dz0 - 2 * SQUISH_CONSTANT_3D;
			}
		}

		/* Contribution (1,1,0) */
		dx3 = dx0 - 1 - 2 * SQUISH_CONSTANT_3D;
		dy3 = dy0 - 1 - 2 * SQUISH_CONSTANT_3D;
		dz3 = dz0 - 0 - 2 * SQUISH_CONSTANT_3D;
		attn3 = 2 - dx3 * dx3 - dy3 * dy3 - dz3 * dz3;
		if(attn3 > 0) {
			attn3 *= attn3;
			value += attn3 * attn3 * extrapolate3(ctx, xsb + 1, ysb + 1, zsb + 0, dx3, dy3, dz3);
		}

		/* Contribution (1,0,1) */
		dx2 = dx3;
		dy2 = dy0 - 0 - 2 * SQUISH_CONSTANT_3D;
		dz2 = dz0 - 1 - 2 * SQUISH_CONSTANT_3D;
		attn2 = 2 - dx2 * dx2 - dy2 * dy2 - dz2 * dz2;
		if(attn2 > 0) {
			attn2 *= attn2;
			value += attn2 * attn2 * extrapolate3(ctx, xsb + 1, ysb + 0, zsb + 1, dx2, dy2, dz2);
		}

		/* Contribution (0,1,1) */
		dx1 = dx0 - 0 - 2 * SQUISH_CONSTANT_3D;
		dy1 = dy3;
		dz1 = dz2;
		attn1 = 2 - dx1 * dx1 - dy1 * dy1 - dz1 * dz1;
		if(attn1 > 0) {
			attn1 *= attn1;
			value += attn1 * attn1 * extrapolate3(ctx, xsb + 0, ysb + 1, zsb + 1, dx1, dy1, dz1);
		}

		/* Contribution (1,1,1) */
		dx0 = dx0 - 1 - 3 * SQUISH_CONSTANT_3D;
		dy0 = dy0 - 1 - 3 * SQUISH_CONSTANT_3D;
		dz0 = dz0 - 1 - 3 * SQUISH_CONSTANT_3D;
		attn0 = 2 - dx0 * dx0 - dy0 * dy0 - dz0 * dz0;
		if(attn0 > 0) {
			attn0 *= attn0;
			value += attn0 * attn0 * extrapolate3(ctx, xsb + 1, ysb + 1, zsb + 1, dx0, dy0, dz0);
		}
	}
	else { /* We're inside the octahedron (Rectified 3-Simplex) in between.
		   Decide between point (0,0,1) and (1,1,0) as closest */
		p1 = xins + yins;
		if(p1 > 1) {
			aScore = p1 - 1;
			aPoint = 0x03;
			aIsFurtherSide = 1;
		}
		else {
			aScore = 1 - p1;
			aPoint = 0x04;
			aIsFurtherSide = 0;
		}

		/* Decide between point (0,1,0) and (1,0,1) as closest */
		p2 = xins + zins;
		if(p2 > 1) {
			bScore = p2 - 1;
			bPoint = 0x05;
			bIsFurtherSide = 1;
		}
		else {
			bScore = 1 - p2;
			bPoint = 0x02;
			bIsFurtherSide = 0;
		}

		/* The closest out of the two (1,0,0) and (0,1,1) will replace the furthest out of the two decided above, if closer. */
		p3 = yins + zins;
		if(p3 > 1) {
			score = p3 - 1;
			if(aScore <= bScore && aScore < score) {
				aScore = score;
				aPoint = 0x06;
				aIsFurtherSide = 1;
			}
			else if(aScore > bScore && bScore < score) {
				bScore = score;
				bPoint = 0x06;
				bIsFurtherSide = 1;
			}
		}
		else {
			score = 1 - p3;
			if(aScore <= bScore && aScore < score) {
				aScore = score;
				aPoint = 0x01;
				aIsFurtherSide = 0;
			}
			else if(aScore > bScore && bScore < score) {
				bScore = score;
				bPoint = 0x01;
				bIsFurtherSide = 0;
			}
		}

		/* Where each of the two closest points are determines how the extra two vertices are calculated. */
		if(aIsFurtherSide == bIsFurtherSide) {
			if(aIsFurtherSide) { /* Both closest points on (1,1,1) side */

				/* One of the two extra points is (1,1,1) */
				dx_ext0 = dx0 - 1 - 3 * SQUISH_CONSTANT_3D;
				dy_ext0 = dy0 - 1 - 3 * SQUISH_CONSTANT_3D;
				dz_ext0 = dz0 - 1 - 3 * SQUISH_CONSTANT_3D;
				xsv_ext0 = xsb + 1;
				ysv_ext0 = ysb + 1;
				zsv_ext0 = zsb + 1;

				/* Other extra point is based on the shared axis. */
				c = (int8_t)(aPoint & bPoint);
				if((c & 0x01) != 0) {
					dx_ext1 = dx0 - 2 - 2 * SQUISH_CONSTANT_3D;
					dy_ext1 = dy0 - 2 * SQUISH_CONSTANT_3D;
					dz_ext1 = dz0 - 2 * SQUISH_CONSTANT_3D;
					xsv_ext1 = xsb + 2;
					ysv_ext1 = ysb;
					zsv_ext1 = zsb;
				}
				else if((c & 0x02) != 0) {
					dx_ext1 = dx0 - 2 * SQUISH_CONSTANT_3D;
					dy_ext1 = dy0 - 2 - 2 * SQUISH_CONSTANT_3D;
					dz_ext1 = dz0 - 2 * SQUISH_CONSTANT_3D;
					xsv_ext1 = xsb;
					ysv_ext1 = ysb + 2;
					zsv_ext1 = zsb;
				}
				else {
					dx_ext1 = dx0 - 2 * SQUISH_CONSTANT_3D;
					dy_ext1 = dy0 - 2 * SQUISH_CONSTANT_3D;
					dz_ext1 = dz0 - 2 - 2 * SQUISH_CONSTANT_3D;
					xsv_ext1 = xsb;
					ysv_ext1 = ysb;
					zsv_ext1 = zsb + 2;
				}
			}
			else { /* Both closest points on (0,0,0) side */

				/* One of the two extra points is (0,0,0) */
				dx_ext0 = dx0;
				dy_ext0 = dy0;
				dz_ext0 = dz0;
				xsv_ext0 = xsb;
				ysv_ext0 = ysb;
				zsv_ext0 = zsb;

				/* Other extra point is based on the omitted axis. */
				c = (int8_t)(aPoint | bPoint);
				if((c & 0x01) == 0) {
					dx_ext1 = dx0 + 1 - SQUISH_CONSTANT_3D;
					dy_ext1 = dy0 - 1 - SQUISH_CONSTANT_3D;
					dz_ext1 = dz0 - 1 - SQUISH_CONSTANT_3D;
					xsv_ext1 = xsb - 1;
					ysv_ext1 = ysb + 1;
					zsv_ext1 = zsb + 1;
				}
				else if((c & 0x02) == 0) {
					dx_ext1 = dx0 - 1 - SQUISH_CONSTANT_3D;
					dy_ext1 = dy0 + 1 - SQUISH_CONSTANT_3D;
					dz_ext1 = dz0 - 1 - SQUISH_CONSTANT_3D;
					xsv_ext1 = xsb + 1;
					ysv_ext1 = ysb - 1;
					zsv_ext1 = zsb + 1;
				}
				else {
					dx_ext1 = dx0 - 1 - SQUISH_CONSTANT_3D;
					dy_ext1 = dy0 - 1 - SQUISH_CONSTANT_3D;
					dz_ext1 = dz0 + 1 - SQUISH_CONSTANT_3D;
					xsv_ext1 = xsb + 1;
					ysv_ext1 = ysb + 1;
					zsv_ext1 = zsb - 1;
				}
			}
		}
		else { /* One point on (0,0,0) side, one point on (1,1,1) side */
			if(aIsFurtherSide) {
				c1 = aPoint;
				c2 = bPoint;
			}
			else {
				c1 = bPoint;
				c2 = aPoint;
			}

			/* One contribution is a permutation of (1,1,-1) */
			if((c1 & 0x01) == 0) {
				dx_ext0 = dx0 + 1 - SQUISH_CONSTANT_3D;
				dy_ext0 = dy0 - 1 - SQUISH_CONSTANT_3D;
				dz_ext0 = dz0 - 1 - SQUISH_CONSTANT_3D;
				xsv_ext0 = xsb - 1;
				ysv_ext0 = ysb + 1;
				zsv_ext0 = zsb + 1;
			}
			else if((c1 & 0x02) == 0) {
				dx_ext0 = dx0 - 1 - SQUISH_CONSTANT_3D;
				dy_ext0 = dy0 + 1 - SQUISH_CONSTANT_3D;
				dz_ext0 = dz0 - 1 - SQUISH_CONSTANT_3D;
				xsv_ext0 = xsb + 1;
				ysv_ext0 = ysb - 1;
				zsv_ext0 = zsb + 1;
			}
			else {
				dx_ext0 = dx0 - 1 - SQUISH_CONSTANT_3D;
				dy_ext0 = dy0 - 1 - SQUISH_CONSTANT_3D;
				dz_ext0 = dz0 + 1 - SQUISH_CONSTANT_3D;
				xsv_ext0 = xsb + 1;
				ysv_ext0 = ysb + 1;
				zsv_ext0 = zsb - 1;
			}

			/* One contribution is a permutation of (0,0,2) */
			dx_ext1 = dx0 - 2 * SQUISH_CONSTANT_3D;
			dy_ext1 = dy0 - 2 * SQUISH_CONSTANT_3D;
			dz_ext1 = dz0 - 2 * SQUISH_CONSTANT_3D;
			xsv_ext1 = xsb;
			ysv_ext1 = ysb;
			zsv_ext1 = zsb;
			if((c2 & 0x01) != 0) {
				dx_ext1 -= 2;
				xsv_ext1 += 2;
			}
			else if((c2 & 0x02) != 0) {
				dy_ext1 -= 2;
				ysv_ext1 += 2;
			}
			else {
				dz_ext1 -= 2;
				zsv_ext1 += 2;
			}
		}

		/* Contribution (1,0,0) */
		dx1 = dx0 - 1 - SQUISH_CONSTANT_3D;
		dy1 = dy0 - 0 - SQUISH_CONSTANT_3D;
		dz1 = dz0 - 0 - SQUISH_CONSTANT_3D;
		attn1 = 2 - dx1 * dx1 - dy1 * dy1 - dz1 * dz1;
		if(attn1 > 0) {
			attn1 *= attn1;
			value += attn1 * attn1 * extrapolate3(ctx, xsb + 1, ysb + 0, zsb + 0, dx1, dy1, dz1);
		}

		/* Contribution (0,1,0) */
		dx2 = dx0 - 0 - SQUISH_CONSTANT_3D;
		dy2 = dy0 - 1 - SQUISH_CONSTANT_3D;
		dz2 = dz1;
		attn2 = 2 - dx2 * dx2 - dy2 * dy2 - dz2 * dz2;
		if(attn2 > 0) {
			attn2 *= attn2;
			value += attn2 * attn2 * extrapolate3(ctx, xsb + 0, ysb + 1, zsb + 0, dx2, dy2, dz2);
		}

		/* Contribution (0,0,1) */
		dx3 = dx2;
		dy3 = dy1;
		dz3 = dz0 - 1 - SQUISH_CONSTANT_3D;
		attn3 = 2 - dx3 * dx3 - dy3 * dy3 - dz3 * dz3;
		if(attn3 > 0) {
			attn3 *= attn3;
			value += attn3 * attn3 * extrapolate3(ctx, xsb + 0, ysb + 0, zsb + 1, dx3, dy3, dz3);
		}

		/* Contribution (1,1,0) */
		dx4 = dx0 - 1 - 2 * SQUISH_CONSTANT_3D;
		dy4 = dy0 - 1 - 2 * SQUISH_CONSTANT_3D;
		dz4 = dz0 - 0 - 2 * SQUISH_CONSTANT_3D;
		attn4 = 2 - dx4 * dx4 - dy4 * dy4 - dz4 * dz4;
		if(attn4 > 0) {
			attn4 *= attn4;
			value += attn4 * attn4 * extrapolate3(ctx, xsb + 1, ysb + 1, zsb + 0, dx4, dy4, dz4);
		}

		/* Contribution (1,0,1) */
		dx5 = dx4;
		dy5 = dy0 - 0 - 2 * SQUISH_CONSTANT_3D;
		dz5 = dz0 - 1 - 2 * SQUISH_CONSTANT_3D;
		attn5 = 2 - dx5 * dx5 - dy5 * dy5 - dz5 * dz5;
		if(attn5 > 0) {
			attn5 *= attn5;
			value += attn5 * attn5 * extrapolate3(ctx, xsb + 1, ysb + 0, zsb + 1, dx5, dy5, dz5);
		}

		/* Contribution (0,1,1) */
		dx6 = dx0 - 0 - 2 * SQUISH_CONSTANT_3D;
		dy6 = dy4;
		dz6 = dz5;
		attn6 = 2 - dx6 * dx6 - dy6 * dy6 - dz6 * dz6;
		if(attn6 > 0) {
			attn6 *= attn6;
			value += attn6 * attn6 * extrapolate3(ctx, xsb + 0, ysb + 1, zsb + 1, dx6, dy6, dz6);
		}
	}

	/* First extra vertex */
	attn_ext0 = 2 - dx_ext0 * dx_ext0 - dy_ext0 * dy_ext0 - dz_ext0 * dz_ext0;
	if(attn_ext0 > 0)
	{
		attn_ext0 *= attn_ext0;
		value += attn_ext0 * attn_ext0 * extrapolate3(ctx, xsv_ext0, ysv_ext0, zsv_ext0, dx_ext0, dy_ext0, dz_ext0);
	}

	/* Second extra vertex */
	attn_ext1 = 2 - dx_ext1 * dx_ext1 - dy_ext1 * dy_ext1 - dz_ext1 * dz_ext1;
	if(attn_ext1 > 0)
	{
		attn_ext1 *= attn_ext1;
		value += attn_ext1 * attn_ext1 * extrapolate3(ctx, xsv_ext1, ysv_ext1, zsv_ext1, dx_ext1, dy_ext1, dz_ext1);
	}

	return value / NORM_CONSTANT_3D;
}

void ccnGenerateOpenSimplex2D(
	ccnNoise *noise,
	ccnNoiseConfiguration *configuration,
	unsigned int scale)
{
	unsigned int i;
	unsigned int size = noise->width * noise->height;
	unsigned int xSteps = (unsigned int)ceil((float)noise->width / scale);
	unsigned int ySteps = (unsigned int)ceil((float)noise->height / scale);
	unsigned int xOffset = 0;
	unsigned int yOffset = 0;
	float multiplier = configuration->range.high - configuration->range.low;
	struct osn_context *ctx;

	// Calculate offset
	ccnPoint offset = (ccnPoint){ configuration->x * xSteps, configuration->y * ySteps };

	// Calculate subsection offsets
	if(noise->width < scale) {
		offset.x = (int)floor(offset.x * ((float)noise->width / scale));
		xOffset = ccnFloorMod(configuration->x, scale / noise->width) * noise->width;
	}

	if(noise->height < scale) {
		offset.y = (int)floor(offset.y * ((float)noise->height / scale));
		yOffset = ccnFloorMod(configuration->y, scale / noise->height) * noise->height;
	}

	open_simplex_noise(configuration->seed, &ctx);

	float min = 1;
	float max = 0;

	for(i = 0; i < size; ++i) {
		unsigned int y = i / noise->width;
		unsigned int x = i - y * noise->width;
		float v = open_simplex_noise2(ctx, (float)x / scale, (float)y / scale);

		if(v > max) max = v;
		if(v < min) min = v;

		ccnStore(noise->values + i, configuration->storeMethod, v * multiplier + configuration->range.low);
	}

	printf("Min %f\nMax %f\n", min, max);

	open_simplex_noise_free(ctx);
}