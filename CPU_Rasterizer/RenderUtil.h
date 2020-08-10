#include "gz.h"
#include "ExtraDefinition.h"
#include <cmath>

#ifndef RENDER_UTIL
#define RENDER_UTIL

namespace RenderUtil {

	float deg2Rad(float degree);
	void matMult(int size, const GzMatrix a, const GzMatrix b, GzMatrix outMat);
	int coordTransform(int size, const GzMatrix transform, GzCoord4 pos);
	void copyMat(int size, const GzMatrix from, GzMatrix to);
	int normalize(GzCoord coord);
	void getBoundingbox(const GzCoord* vertList, float& minX, float& maxX, float& minY, float& maxY);
	void interpolate3d(float a, float b, float c, float d, float x, float y, float& z);
	int roundToInt(float v);
	void ySort(GzCoord* vertList, int* indexList);
	int LEE(const GzCoord head, const GzCoord tail, float x, float y);
	// v1 -> v2 v2 -> v3 v3 -> v1
	int boundByTri(const GzCoord v1, const GzCoord v2, const GzCoord v3, float x, float y);
	int crossproduct(const GzCoord v1, const GzCoord v2, GzCoord outNormal);
	float dotproduct(const GzCoord v1, const GzCoord v2);
	void IntensityBoundCheck(GzIntensity &inVal);
	void DepthBoundCheck(GzDepth &inVal);
	void reflect(const GzCoord normal, const GzCoord in, GzCoord out);
	void getUrinaryRotation(const GzMatrix in, GzMatrix result);
	int getPlanecoefficient(GzCoord* vertList,
		float &a, float &b, float &c, float &d);
	void getVirtualSpacePoint(const GzCoord* vertList, const GzCoord valueList, GzCoord* resultList);
	void color_x_color(const GzColor c1, const GzColor c2, GzColor result);
	void color_x_color(const GzColor c1, const GzColor c2, const GzColor c3, GzColor result);
	void color_x_float(const GzColor c1, float f, GzColor result);
	void color_plus_color(const GzColor c1, const GzColor c2, GzColor result);
	void color_plus_color(const GzColor c1, const GzColor c2, const GzColor c3, GzColor result);
	//void clamp01(GzColor c) {
	//	if (c[RED] < 0.0f) c[RED] = 0.0f;
	//	if (c[RED] > 1.0f) c[RED] = 1.0f;

	//	if (c[GREEN] < 0.0f) c[GREEN] = 0.0f;
	//	if (c[GREEN] > 1.0f) c[GREEN] = 1.0f;

	//	if (c[BLUE] < 0.0f) c[BLUE] = 0.0f;
	//	if (c[BLUE] > 1.0f) c[BLUE] = 1.0f;
	//}
}
#endif // !RENDER_UTIL
