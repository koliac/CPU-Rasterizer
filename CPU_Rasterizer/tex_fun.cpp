/* Texture functions for cs580 GzLib	*/
#include    "stdafx.h" 
#include	"stdio.h"
#include	"Gz.h"
#include <cmath>

GzColor	*image=NULL;
int xs, ys;
int reset = 1;

int getIndex(int x, int y, int xres) { return x + y * xres; }
/* Image texture function */
int tex_fun(float u, float v, GzColor color)
{
  unsigned char		pixel[3];
  unsigned char     dummy;
  char  		foo[8];
  int   		i, j;
  FILE			*fd;

  if (reset) {          /* open and load texture file */
    fd = fopen ("texture", "rb");
    if (fd == NULL) {
      fprintf (stderr, "texture file not found\n");
      exit(-1);
    }
    fscanf (fd, "%s %d %d %c", foo, &xs, &ys, &dummy);
    image = (GzColor*)malloc(sizeof(GzColor)*(xs+1)*(ys+1));
    if (image == NULL) {
      fprintf (stderr, "malloc for texture image failed\n");
      exit(-1);
    }

    for (i = 0; i < xs*ys; i++) {	/* create array of GzColor values */
      fread(pixel, sizeof(pixel), 1, fd);
      image[i][RED] = (float)((int)pixel[RED]) * (1.0 / 255.0);
      image[i][GREEN] = (float)((int)pixel[GREEN]) * (1.0 / 255.0);
      image[i][BLUE] = (float)((int)pixel[BLUE]) * (1.0 / 255.0);
      }

    reset = 0;          /* init is done */
	fclose(fd);
  }

/* bounds-test u,v to make sure nothing will overflow image array bounds */
/* determine texture cell corner values and perform bilinear interpolation */
/* set color to interpolated GzColor value and return */

  // bound check - clip the texture [0, 1];
  /*float desiredU = u ;
  if (desiredU > 1.0f) desiredU = 1.0f;
  if (desiredU < 0.0f) desiredU = 0.0f;
  float desiredV = v;
  if (desiredV > 1.0f) desiredV = 1.0f;
  if (desiredV < 0.0f) desiredV = 0.0f;*/

  // apply bilinear filter
  float s = u*(static_cast<float>(xs)-1.0f);
  float t = v*(static_cast<float>(ys)-1.0f);

  if (s > (float)xs - 1.0f) s = (float)xs - 1.0f;
  if (s < 0.0) s = 0.0f;

  if (t > (float)ys - 1.0f) t = (float)ys - 1.0f;
  if (t < 0.0) t = 0.0f;

  int floorS = (int)std::floor(s);
  int ceilS = (int)std::ceil(s);

  int floorT = (int)std::floor(t);
  int ceilT = (int)std::ceil(t);

  GzColor A{ 
	  image[getIndex(floorS, floorT, xs)][RED],
	  image[getIndex(floorS, floorT, xs)][GREEN],
	  image[getIndex(floorS, floorT, xs)][BLUE] 
  };

  GzColor B{
	 image[getIndex(ceilS, floorT, xs)][RED],
	 image[getIndex(ceilS, floorT, xs)][GREEN],
	 image[getIndex(ceilS, floorT, xs)][BLUE]
  };

  GzColor C{
	 image[getIndex(ceilS, ceilT, xs)][RED],
	 image[getIndex(ceilS, ceilT, xs)][GREEN],
	 image[getIndex(ceilS, ceilT, xs)][BLUE],
  };

  GzColor D{
	 image[getIndex(floorS, ceilT, xs)][RED],
	 image[getIndex(floorS, ceilT, xs)][GREEN],
	 image[getIndex(floorS, ceilT, xs)][BLUE]
  };

  float ds = s - std::floor(s);
  float dt = t - std::floor(t);

  color[RED] = (ds * dt)*C[RED] + ((1 - ds)*dt)*D[RED] + (ds * (1 - dt))*B[RED] + ((1 - ds)*(1 - dt))*A[RED];
  color[GREEN] = (ds * dt)*C[GREEN] + ((1 - ds)*dt)*D[GREEN] + (ds * (1 - dt))*B[GREEN] + ((1 - ds)*(1 - dt))*A[GREEN];
  color[BLUE] = (ds * dt)*C[BLUE] + ((1 - ds)*dt)*D[BLUE] + (ds * (1 - dt))*B[BLUE] + ((1 - ds)*(1 - dt))*A[BLUE];

  //color[RED] = A[RED];
  //color[GREEN] = A[GREEN];
  //color[BLUE] = A[BLUE];

  //color[RED] = desiredU;
  //color[GREEN] = desiredV;
  //color[BLUE] = 0.0f;

  return GZ_SUCCESS;
}
void voronoiNoise(float u, float v, GzColor result) {

	float baseU = std::floor(u);
	float baseV = std::floor(v);
	float du = 0.5f*std::sin(u) + 0.5f;
	float dv = 0.5f*std::sin(v) + 0.5f;
	float mindist = 10.0f;
	float dist;
	

	for (int x = -1; x <= 1; ++x)
	{
		for (int y = -1;y <= 1; ++y)
		{
			float cellU = baseU + x;
			float cellV = baseV + y;
			float cellPosU = cellU+du;
			float cellPosV = cellV+dv;
			float toCellU = cellPosU - u;
			float toCellV = cellPosV - v;
			dist = std::sqrt(toCellU*toCellU + toCellV * toCellV);
			if (dist < mindist) {
				mindist = dist;
			}
		}
	}
	result[RED] = mindist;
	result[GREEN] = mindist;
	result[BLUE] = mindist;
}
/* Procedural texture function */
int ptex_fun(float u, float v, GzColor color)
{

	float desiredU = u;
	if (desiredU > 1.0f) desiredU = 1.0f;
	if (desiredU < 0.0f) desiredU = 0.0f;
	float desiredV = v;
	if (desiredV > 1.0f) desiredV = 1.0f;
	if (desiredV < 0.0f) desiredV = 0.0f;
	voronoiNoise(desiredU/0.1f,desiredV/0.1f, color);

	
	return GZ_SUCCESS;
}

/* Free texture memory */
int GzFreeTexture()
{
	if(image!=NULL)
		free(image);
	return GZ_SUCCESS;
}

