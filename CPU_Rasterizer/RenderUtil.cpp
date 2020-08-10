#include "stdafx.h"
#include "RenderUtil.h"
#include "ExtraDefinition.h"
#include <cmath>
#include <utility>
#include "rend.h"

namespace RenderUtil {

	float deg2Rad(float degree) {
		return (degree / 180.0f)*PI;
	}

	void matMult(int size, const GzMatrix a, const GzMatrix b, GzMatrix outMat) {

		for (unsigned int i = 0; i < size; ++i)
		{
			for (unsigned int j = 0; j < size; ++j)
			{
				outMat[i][j] = 0.0f;

			}
		}
		for (unsigned int i = 0; i < size; ++i)
		{
			for (unsigned int j = 0; j < size; ++j)
			{
				for (unsigned int k = 0; k < size; ++k)
				{
					outMat[i][j] += a[i][k] * b[k][j];
				}
			}
		}

	}

	int coordTransform(int size, const GzMatrix transform, GzCoord4 pos) {
		if (size != 4) return GZ_FAILURE;
		GzCoord4 result{ 0.0f,0.0f,0.0f,0.0f };

		for (unsigned int i = 0; i < size; ++i)
		{
			for (unsigned int j = 0; j < size; ++j)
			{
				result[i] += transform[i][j] * pos[j];
			}
		}

		pos[X] = result[X];
		pos[Y] = result[Y];
		pos[Z] = result[Z];
		pos[W] = result[W];

		return GZ_SUCCESS;
	}

	void copyMat(int size, const GzMatrix from, GzMatrix to) {
		for (unsigned int i = 0; i < size; ++i)
		{
			for (unsigned int j = 0; j < size; ++j)
			{
				to[i][j] = from[i][j];
			}
		}
	}

	int normalize(GzCoord coord) {
		float u = std::sqrt(coord[X] * coord[X] + coord[Y] * coord[Y] + coord[Z] * coord[Z]);
		if (u == 0.0f) {
			return GZ_FAILURE;
		}
		coord[X] = coord[X] / u;
		coord[Y] = coord[Y] / u;
		coord[Z] = coord[Z] / u;
		return GZ_SUCCESS;
	}

	void getBoundingbox(const GzCoord* vertList, float& minX, float& maxX, float& minY, float& maxY) {
		minX = vertList[0][X];
		maxX = vertList[0][X];
		minY = vertList[0][Y];
		maxY = vertList[0][Y];

		if (vertList[1][X] < minX) minX = vertList[1][X];
		if (vertList[1][X] > maxX) maxX = vertList[1][X];
		if (vertList[1][Y] < minY) minY = vertList[1][Y];
		if (vertList[1][Y] > maxY) maxY = vertList[1][Y];

		if (vertList[2][X] < minX) minX = vertList[2][X];
		if (vertList[2][X] > maxX) maxX = vertList[2][X];
		if (vertList[2][Y] < minY) minY = vertList[2][Y];
		if (vertList[2][Y] > maxY) maxY = vertList[2][Y];


	}

	void interpolate3d(float a, float b, float c, float d, float x, float y, float& z) {
		z = ((a*x + b * y + d) / (-c));
	}

	int roundToInt(float v) {
		return std::floor(v);
	}

	void ySort(GzCoord* vertList, int* indexList) {

		int rtn[3]{ 0,1,2 };
		if (vertList[0][Y] > vertList[1][Y]) {
			std::swap(vertList[0], vertList[1]);
			std::swap(rtn[0], rtn[1]);

		}
		if (vertList[0][Y] > vertList[2][Y]) {
			std::swap(vertList[0], vertList[2]);
			std::swap(rtn[0], rtn[2]);

		}
		if (vertList[1][Y] > vertList[2][Y]) {
			std::swap(vertList[1], vertList[2]);
			std::swap(rtn[1], rtn[2]);

		}
		indexList[0] = rtn[0];
		indexList[1] = rtn[1];
		indexList[2] = rtn[2];


	}

	int LEE(const GzCoord head, const GzCoord tail, float x, float y) {
		float dx{ head[X] - tail[X] };
		float dy{ head[Y] - tail[Y] };

		float a{ dy };
		float b{ -dx };
		float c{ (dx*tail[Y]) - (dy*tail[X]) };
		float result{ (a*x) + (b*y) + c };
		if (result < 0.0f) {
			return R;
		}
		else {
			return L;
		}
	}

	int boundByTri(const GzCoord v1, const GzCoord v2, const GzCoord v3, float x, float y) {
		if ((LEE(v2, v1, x, y) == LEE(v3, v2, x, y))
			&& (LEE(v3, v2, x, y) == LEE(v1, v3, x, y))) {
			return 1;
		}
		else {
			return 0;
		}
	}

	int crossproduct(const GzCoord v1, const GzCoord v2, GzCoord outNormal) {
		int status{ GZ_SUCCESS };
		outNormal[0] = v1[1] * v2[2] - v1[2] * v2[1];
		outNormal[1] = v1[2] * v2[0] - v1[0] * v2[2];
		outNormal[2] = v1[0] * v2[1] - v1[1] * v2[0];
		return status;
	}

	float dotproduct(const GzCoord v1, const GzCoord v2) {
		return v1[X] * v2[X] + v1[Y] * v2[Y] + v1[Z] * v2[Z];
	}
	void IntensityBoundCheck(GzIntensity &inVal) {
		if (inVal < 0) inVal = 0;
		if (inVal > MAXINTENSITY) inVal = MAXINTENSITY;
	}

	void DepthBoundCheck(GzDepth &inVal) {
		if (inVal < 0) inVal = 0;
		if (inVal > INT_MAX) inVal = INT_MAX;
	}
	void reflect(const GzCoord normal, const GzCoord in, GzCoord out) {
		float ndotl = dotproduct(normal, in);
		out[X] = 2.0f * ndotl*normal[X] - in[X];
		out[Y] = 2.0f *ndotl*normal[Y] - in[Y];
		out[Z] = 2.0f*ndotl*normal[Z] - in[Z];
	}

	void getUrinaryRotation(const GzMatrix in, GzMatrix result) {
		float kx = std::sqrt(
			in[0][0] * in[0][0] +
			in[0][1] * in[0][1] +
			in[0][2] * in[0][2]);
		float ky = std::sqrt(
			in[1][0] * in[1][0] +
			in[1][1] * in[1][1] +
			in[1][2] * in[1][2]);
		float kz = std::sqrt(
			in[2][0] * in[2][0] +
			in[2][1] * in[2][1] +
			in[2][2] * in[2][2]);

		//x
		result[0][0] = in[0][0] / kx;
		result[0][1] = in[0][1] / kx;
		result[0][2] = in[0][2] / kx;

		//y
		result[1][0] = in[1][0] / ky;
		result[1][1] = in[1][1] / ky;
		result[1][2] = in[1][2] / ky;

		//z
		result[2][0] = in[2][0] / kz;
		result[2][1] = in[2][1] / kz;
		result[2][2] = in[2][2] / kz;

		result[3][0] = 0.0f;
		result[3][1] = 0.0f;
		result[3][2] = 0.0f;
		result[3][3] = 1.0f;

		result[0][3] = 0.0f;
		result[1][3] = 0.0f;
		result[2][3] = 0.0f;

	}
	int getPlanecoefficient(GzCoord* vertList,
		float &a, float &b, float &c, float &d) {
		//RenderUtil::ySort(vertList);
		int test{ RenderUtil::LEE(vertList[0],vertList[2],vertList[1][X],vertList[1][Y]) };
		GzCoord norm;
		int status;
		if (test == L) {
			status = RenderUtil::crossproduct(
				GzCoord{ vertList[1][X] - vertList[0][X],
				vertList[1][Y] - vertList[0][Y],
				vertList[1][Z] - vertList[0][Z] },
				GzCoord{ vertList[2][X] - vertList[0][X],
				vertList[2][Y] - vertList[0][Y],
				vertList[2][Z] - vertList[0][Z] },
				norm);
		}
		else {
			status = RenderUtil::crossproduct(
				GzCoord{ vertList[2][X] - vertList[0][X],
				vertList[2][Y] - vertList[0][Y],
				vertList[2][Z] - vertList[0][Z] },
				GzCoord{ vertList[1][X] - vertList[0][X],
				vertList[1][Y] - vertList[0][Y],
				vertList[1][Z] - vertList[0][Z] },
				norm);
		}
		//clockwise


		if (status == GZ_FAILURE) {
			return GZ_FAILURE;
		}

		a= norm[0];
		b= norm[1];
		c= norm[2];
		double ax{ a*vertList[1][X] };
		double by{ b*vertList[1][Y] };
		double cz{ c*vertList[1][Z] };
		d= static_cast<float>(std::round(-(ax + by + cz)));
		return GZ_SUCCESS;
	}

	void getVirtualSpacePoint(const GzCoord* vertList, const GzCoord valueList, GzCoord* resultList) {
		resultList[0][X] = vertList[0][X];
		resultList[0][Y] = vertList[0][Y];
		resultList[0][Z] = valueList[0];

		resultList[1][X] = vertList[1][X];
		resultList[1][Y] = vertList[1][Y];
		resultList[1][Z] = valueList[1];

		resultList[2][X] = vertList[2][X];
		resultList[2][Y] = vertList[2][Y];
		resultList[2][Z] = valueList[2];

	}

	void color_x_color(const GzColor c1, const GzColor c2, GzColor result)
	{
		result[RED] = c1[RED] * c2[RED];
		result[GREEN] = c1[GREEN] * c2[GREEN];
		result[BLUE] = c1[BLUE] * c2[BLUE];
	}

	void color_x_color(const GzColor c1, const GzColor c2, const GzColor c3, GzColor result)
	{
		result[RED] = c1[RED] * c2[RED]*c3[RED];
		result[GREEN] = c1[GREEN] * c2[GREEN]*c3[GREEN];
		result[BLUE] = c1[BLUE] * c2[BLUE]*c3[BLUE];
	}

	void color_x_float(const GzColor c1, float f, GzColor result)
	{
		result[RED] = c1[RED] * f;
		result[GREEN] = c1[GREEN] * f;
		result[BLUE] = c1[BLUE] * f;
	}

	void color_plus_color(const GzColor c1, const GzColor c2, GzColor result)
	{
		result[RED] = c1[RED] + c2[RED];
		result[GREEN] = c1[GREEN] + c2[GREEN];
		result[BLUE] = c1[BLUE] + c2[BLUE];
	}

	void color_plus_color(const GzColor c1, const GzColor c2, const GzColor c3, GzColor result)
	{
		result[RED] = c1[RED] + c2[RED] + c3[RED];
		result[GREEN] = c1[GREEN] + c2[GREEN] + c3[GREEN];
		result[BLUE] = c1[BLUE] + c2[BLUE] + c3[BLUE];
	}


}
