/* CS580 Homework 3 */

#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"
#include <utility>
#include <string>
#include "ExtraDefinition.h"
#include "RenderUtil.h"

#ifndef PI
#define PI (float) 3.14159265358979323846
#endif // PI
int shadeSurface(const GzColor Ka, const GzColor Ks, const GzColor Kd, const float spec,
	const GzCoord normal, GzCoord eye, GzLight* lights, const GzLight ambient, GzColor outColor) {
	bool status = GZ_SUCCESS;
	GzColor specular{ 0.0f,0.0f,0.0f };
	GzColor diffuse{ 0.0f,0.0f,0.0f };
	GzColor ambientFinal;
	RenderUtil::color_x_color(Ka, ambient.color, ambientFinal);
	for (unsigned int li = 0; li < 3; ++li)
	{



		status |= RenderUtil::normalize(lights[li].direction);
		status |= RenderUtil::normalize(eye);
		//RenderUtil::reflect(normalListImg[vi], this->lights[li].direction, r);
		// eye vector 


		float ndote = RenderUtil::dotproduct(normal, eye);
		float ndotl = RenderUtil::dotproduct(normal, lights[li].direction);

		GzCoord r;
		if (ndotl > 0.0f && ndote > 0.0f) {
			RenderUtil::reflect(normal, lights[li].direction, r);
		}
		else if (ndotl < 0.0f &&ndote < 0.0f) {
			GzCoord newN{ -normal[X], -normal[Y], -normal[Z] };
			ndotl = RenderUtil::dotproduct(newN, lights[li].direction);
			RenderUtil::reflect(newN, lights[li].direction, r);
		}
		else {
			continue;
		}

		float rdote = RenderUtil::dotproduct(r, eye);
		if (rdote > 1.0f) rdote = 1.0f;
		if (rdote < 0.0f) rdote = 0.0f;
		rdote = std::pow(rdote, spec);

		specular[RED] += lights[li].color[RED] * rdote;
		specular[GREEN] += lights[li].color[GREEN] * rdote;
		specular[BLUE] += lights[li].color[BLUE] * rdote;

		diffuse[RED] += lights[li].color[RED] * ndotl;
		diffuse[GREEN] += lights[li].color[GREEN] * ndotl;
		diffuse[BLUE] += lights[li].color[BLUE] * ndotl;
	}

	GzColor specularFinal;
	GzColor diffuseFinal;
	RenderUtil::color_x_color(Ks, specular, specularFinal);
	RenderUtil::color_x_color(Kd, diffuse, diffuseFinal);
	RenderUtil::color_plus_color(specularFinal, diffuseFinal, ambientFinal, outColor);
	return status;
}
int GzRender::GzRotXMat(float degree, GzMatrix mat)
{
/* HW 3.1
// Create rotate matrix : rotate along x axis
// Pass back the matrix using mat value
*/
	float rad = RenderUtil::deg2Rad(degree);
	GzMatrix rotX{
		{1.0f,0.0f,0.0f,0.0f},
		{0.0f,std::cos(rad),-std::sin(rad), 0.0f},
		{0.0f,std::sin(rad),std::cos(rad), 0.0f},
		{0.0f,0.0f,0.0f,1.0f}
	};
	RenderUtil::copyMat(4, rotX, mat);
	return GZ_SUCCESS;
}

int GzRender::GzRotYMat(float degree, GzMatrix mat)
{
/* HW 3.2
// Create rotate matrix : rotate along y axis
// Pass back the matrix using mat value
*/
	float rad = RenderUtil::deg2Rad(degree);
	GzMatrix rotY{
		{std::cos(rad),0.0f,std::sin(rad), 0.0f},
		{0.0f,1.0f,0.0f,0.0f},
		{-std::sin(rad),0.0f,std::cos(rad), 0.0f},
		{0.0f,0.0f,0.0f,1.0f}
	};
	RenderUtil::copyMat(4, rotY, mat);
	return GZ_SUCCESS;
}

int GzRender::GzRotZMat(float degree, GzMatrix mat)
{
/* HW 3.3
// Create rotate matrix : rotate along z axis
// Pass back the matrix using mat value
*/
	float rad = RenderUtil::deg2Rad(degree);
	GzMatrix rotZ{
		{std::cos(rad),-std::sin(rad),0.0f,0.0f},
		{std::sin(rad),std::cos(rad),0.0f,0.0f},
		{0.0f,0.0f,1.0f,0.0f},
		{0.0f,0.0f,0.0f,1.0f}
	};
	RenderUtil::copyMat(4, rotZ, mat);
	return GZ_SUCCESS;
}

int GzRender::GzTrxMat(GzCoord translate, GzMatrix mat)
{
/* HW 3.4
// Create translation matrix
// Pass back the matrix using mat value
*/
	GzMatrix trx{
			{1.0f,0.0f,0.0f,translate[X]},
			{0.0f,1.0f,0.0f,translate[Y]},
			{0.0f,0.0f,1.0f,translate[Z]},
			{0.0f,0.0f,0.0f,1.0f}
	};

	RenderUtil::copyMat(4, trx, mat);
	return GZ_SUCCESS;
}


int GzRender::GzScaleMat(GzCoord scale, GzMatrix mat)
{
/* HW 3.5
// Create scaling matrix
// Pass back the matrix using mat value
*/
	GzMatrix scaleMat{
		{scale[X],0.0f,0.0f,0.0f},
		{0.0f,scale[Y],0.0f,0.0f},
		{0.0f,0.0f,scale[Z],0.0f},
		{0.0f,0.0f,0.0f,1.0f}
	};

	RenderUtil::copyMat(4, scaleMat, mat);
	return GZ_SUCCESS;
}


GzRender::GzRender(int xRes, int yRes, float xshiftScalar, float yshiftScalar)
{
/* HW1.1 create a framebuffer for MS Windows display:
 -- set display resolution
 -- allocate memory for framebuffer : 3 bytes(b, g, r) x width x height
 -- allocate memory for pixel buffer
 */
	int desiredXRes{ xRes > MAXXRES ? MAXXRES : xRes };
	int desiredYRes{ yRes > MAXYRES ? MAXYRES : yRes };

	xres = static_cast<unsigned short>(desiredXRes);
	yres = static_cast<unsigned short>(desiredYRes);
	this->xshiftScalar = xshiftScalar;
	this->yshiftScalar = yshiftScalar;
	matlevel = -1;
	normMatlevel = -1;

	pixelbuffer = new GzPixel[xres*yres];
	framebuffer = (char*)malloc(3 * sizeof(char) * xRes * yRes);

	this->m_camera.position[X] = DEFAULT_IM_X;
	this->m_camera.position[Y] = DEFAULT_IM_Y;
	this->m_camera.position[Z] = DEFAULT_IM_Z;

	this->m_camera.lookat[X] = 0.0f;
	this->m_camera.lookat[Y] = 0.0f;
	this->m_camera.lookat[Z] = 0.0f;

	this->m_camera.worldup[X] = 0.0f;
	this->m_camera.worldup[Y] = 1.0f;
	this->m_camera.worldup[Z] = 0.0f;

	this->m_camera.FOV = DEFAULT_FOV;

	GzMatrix p2s{
		{xres / 2.0f,0.0f,0.0f,xres / 2.0f},
		{0.0f,-yres / 2.0f,0.0f,yres / 2.0f},
		{0.0f,0.0f,INT_MAX,0.0f},
		{0.0f,0.0f,0.0f,1.0f}
	};
	RenderUtil::copyMat(4, p2s, this->Xsp);

/* HW 3.6
- setup Xsp and anything only done once 
- init default camera 
*/ 
}

GzRender::~GzRender()
{
/* HW1.2 clean up, free buffer memory */
	delete[] pixelbuffer;
	free(framebuffer);

	pixelbuffer = NULL;
	framebuffer = NULL;
}

int GzRender::GzDefault()
{
/* HW1.3 set pixel buffer to some default values - start a new frame */

	for (unsigned int i = 0; i < xres*yres; i++)
	{
		// 255, 248, 171
		GzIntensity r{ 255 << 4 };
		GzIntensity g{ 248 << 4 };
		GzIntensity b{ 171 << 4 };
		GzIntensity a{ 1 };
		GzDepth z{ INT_MAX };

		pixelbuffer[i] = GzPixel{ r,g,b,a,z };
	}

	return GZ_SUCCESS;
}

int GzRender::GzBeginRender()
{
/* HW 3.7 
- setup for start of each frame - init frame buffer color,alpha,z
- compute Xiw and projection xform Xpi from camera definition 
- init Ximage - put Xsp at base of stack, push on Xpi and Xiw 
- now stack contains Xsw and app can push model Xforms when needed 
*/ 
	int rtnStatus = GZ_SUCCESS;
	GzDefault();
	matlevel = -1;
	normMatlevel = -1;
	GzCoord forward{
	 this->m_camera.lookat[X] - this->m_camera.position[X],
	 this->m_camera.lookat[Y] - this->m_camera.position[Y],
	 this->m_camera.lookat[Z] - this->m_camera.position[Z]
	};

	rtnStatus |= RenderUtil::normalize(forward);

	float d = RenderUtil::dotproduct(this->m_camera.worldup, forward);
	GzCoord up{ this->m_camera.worldup[X] - (forward[X] * d),
	this->m_camera.worldup[Y] - (forward[Y] * d),
	this->m_camera.worldup[Z] - (forward[Z] * d) };

	rtnStatus |= RenderUtil::normalize(up);

	GzCoord right;
	RenderUtil::crossproduct(up, forward, right);
	float xc = RenderUtil::dotproduct(right, this->m_camera.position);
	float yc = RenderUtil::dotproduct(up, this->m_camera.position);
	float zc = RenderUtil::dotproduct(forward, this->m_camera.position);

	GzMatrix Xiw{
		{right[X],right[Y],right[Z],-xc},
		{up[X],up[Y],up[Z],-yc},
		{forward[X],forward[Y],forward[Z],-zc},
		{0.0f,0.0f,0.0f,1.0f}
	};

	//matMult(4, coordSys, camTraslate, Xiw);

	RenderUtil::copyMat(4, Xiw, this->m_camera.Xiw);

	float denom = std::tan(RenderUtil::deg2Rad(this->m_camera.FOV) / 2.0f);

	GzMatrix Xpi{
		{1.0f,0.0f,0.0f,0.0f},
		{0.0f,1.0f,0.0f,0.0f},
		{0.0f,0.0f,denom,0.0f},
		{0.0f,0.0f,denom,1.0f}

	};

	RenderUtil::copyMat(4, Xpi, this->m_camera.Xpi);


	return GZ_SUCCESS;
}

int GzRender::GzPutCamera(GzCamera camera)
{
/* HW 3.8 
/*- overwrite renderer camera structure with new camera definition
*/
	this->m_camera = camera;
	return GZ_SUCCESS;
}

int GzRender::GzPushMatrix(GzMatrix	matrix)
{
/* HW 3.9 
- push a matrix onto the Ximage stack
- check for stack overflow
*/
	if (matlevel + 1 >= MATLEVELS) {
		return GZ_FAILURE;
	}
	if (matlevel + 1 == 0) {
		RenderUtil::copyMat(4, matrix, Ximage[matlevel + 1]);
	}
	else {
		RenderUtil::matMult(4, Ximage[matlevel], matrix, Ximage[matlevel + 1]);
	}

	++matlevel;
	return GZ_SUCCESS;
}

int GzRender::GzPopMatrix()
{
/* HW 3.10
- pop a matrix off the Ximage stack
- check for stack underflow
*/
	if (matlevel - 1 < 0) {
		return GZ_FAILURE;
	}
	--matlevel;
	return GZ_SUCCESS;
}

int GzRender::GzPushNormalMatrix(GzMatrix matrix) {
	if (normMatlevel + 1 >= MATLEVELS) {
		return GZ_FAILURE;
	}
	if (normMatlevel + 1 == 0) {
		RenderUtil::copyMat(4, matrix, Xnorm[normMatlevel + 1]);
	}
	else {
		RenderUtil::matMult(4, Xnorm[normMatlevel], matrix, Xnorm[normMatlevel + 1]);
	}

	++normMatlevel;
	return GZ_SUCCESS;
}
int GzRender::GzPopNormalMatrix(GzMatrix matrix) {
	if (normMatlevel - 1 < 0) {
		return GZ_FAILURE;
	}
	--normMatlevel;
	return GZ_SUCCESS;
}

int GzRender::GzPut(int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z)
{
/* HW1.4 write pixel values into the buffer */
	if (i < 0 ||
		i > xres - 1 ||
		j < 0 ||
		j > yres - 1) {
		return GZ_FAILURE;
	}


	RenderUtil::IntensityBoundCheck(r);
	RenderUtil::IntensityBoundCheck(g);
	RenderUtil::IntensityBoundCheck(b);
	RenderUtil::IntensityBoundCheck(a);
	RenderUtil::DepthBoundCheck(z);

	GzPixel pixelToSet{ r,g,b,a,z };
	pixelbuffer[ARRAY(i, j)] = pixelToSet;
	return GZ_SUCCESS;
}


int GzRender::GzGet(int i, int j, GzIntensity *r, GzIntensity *g, GzIntensity *b, GzIntensity *a, GzDepth *z)
{
/* HW1.5 retrieve a pixel information from the pixel buffer */
	if (i < 0 ||
		i > xres - 1 ||
		j < 0 ||
		j > yres - 1) {
		return GZ_FAILURE;
	}

	GzPixel p{ pixelbuffer[ARRAY(i, j)] };
	*r = p.red;
	*g = p.green;
	*b = p.blue;
	*a = p.alpha;
	*z = p.z;
	return GZ_SUCCESS;
}


int GzRender::GzFlushDisplay2File(FILE* outfile)
{
/* HW1.6 write image to ppm file -- "P6 %d %d 255\r" */
	if (outfile == nullptr) return GZ_FAILURE;
	std::string header{ "P6" };
	header += " " + std::to_string(xres);
	header += " " + std::to_string(yres);
	header += " 255\r";
	fputs(header.c_str(), outfile);
	for (unsigned int i = 0; i < xres*yres; i++)
	{
		fputc(static_cast<char>(pixelbuffer[i].red >> 4), outfile);
		fputc(static_cast<char>(pixelbuffer[i].green >> 4), outfile);
		fputc(static_cast<char>(pixelbuffer[i].blue >> 4), outfile);
	}
	return GZ_SUCCESS;
}

int GzRender::GzFlushDisplay2FrameBuffer()
{
/* HW1.7 write pixels to framebuffer: 
	- put the pixels into the frame buffer
	- CAUTION: when storing the pixels into the frame buffer, the order is blue, green, and red 
	- NOT red, green, and blue !!!
*/
	unsigned int j{ 0 };
	for (unsigned int i = 0; i < xres*yres; i++)
	{
		framebuffer[j] = static_cast<char>(pixelbuffer[i].blue >> 4);
		framebuffer[j + 1] = static_cast<char>(pixelbuffer[i].green >> 4);
		framebuffer[j + 2] = static_cast<char>(pixelbuffer[i].red >> 4);
		j += 3;
	}
	return GZ_SUCCESS;
}


/***********************************************/
/* HW2 methods: implement from here */

int GzRender::GzPutAttribute(int numAttributes, GzToken	*nameList, GzPointer *valueList) 
{
/* HW 2.1
-- Set renderer attribute states (e.g.: GZ_RGB_COLOR default color)
-- In later homeworks set shaders, interpolaters, texture maps, and lights
*/

/*
- GzPutAttribute() must accept the following tokens/values:

- GZ_RGB_COLOR					GzColor		default flat-shade color
- GZ_INTERPOLATE				int			shader interpolation mode
- GZ_DIRECTIONAL_LIGHT			GzLight
- GZ_AMBIENT_LIGHT            	GzLight		(ignore direction)
- GZ_AMBIENT_COEFFICIENT		GzColor		Ka reflectance
- GZ_DIFFUSE_COEFFICIENT		GzColor		Kd reflectance
- GZ_SPECULAR_COEFFICIENT       GzColor		Ks coef's
- GZ_DISTRIBUTION_COEFFICIENT   float		spec power
*/
	for (unsigned int i = 0; i < numAttributes; ++i) {
		GzToken token = nameList[i];
		switch (token)
		{
		case GZ_RGB_COLOR: {
			GzColor* colPtr{ static_cast<GzColor*>(valueList[i]) };
			float r{ (*colPtr)[RED] };
			float g{ (*colPtr)[GREEN] };
			float b{ (*colPtr)[BLUE] };
			flatcolor[0] = r;
			flatcolor[1] = g;
			flatcolor[2] = b;
			break;
		}
		case GZ_INTERPOLATE: {
			this->interp_mode = *static_cast<int*>(valueList[i]);
			break;
		}
		case GZ_DIRECTIONAL_LIGHT: {
			GzLight* lightsPtr{ static_cast<GzLight*>(valueList[i]) };
			this->lights[i] = *lightsPtr;
			break;
		}
		case GZ_AMBIENT_LIGHT: {
			this->ambientlight = *static_cast<GzLight*>(valueList[i]);
			break;
		}
		case GZ_AMBIENT_COEFFICIENT: {
			GzColor* colPtr{ static_cast<GzColor*>(valueList[i]) };
			this->Ka[RED] = (*colPtr)[RED];
			this->Ka[GREEN] = (*colPtr)[GREEN];
			this->Ka[BLUE] = (*colPtr)[BLUE];
			break;
		}
		case GZ_DIFFUSE_COEFFICIENT: {
			GzColor* colPtr{ static_cast<GzColor*>(valueList[i]) };
			this->Kd[RED] = (*colPtr)[RED];
			this->Kd[GREEN] = (*colPtr)[GREEN];
			this->Kd[BLUE] = (*colPtr)[BLUE];
			break;
		}
		case GZ_SPECULAR_COEFFICIENT: {
			GzColor* colPtr{ static_cast<GzColor*>(valueList[i]) };
			this->Ks[RED] = (*colPtr)[RED];
			this->Ks[GREEN] = (*colPtr)[GREEN];
			this->Ks[BLUE] = (*colPtr)[BLUE];
			break;
		}
		case GZ_DISTRIBUTION_COEFFICIENT: {
			this->spec = *static_cast<float*>(valueList[i]);
			break;
		}
		case GZ_TEXTURE_MAP: {
			this->tex_fun= static_cast<GzTexture>(valueList[i]);
			break;
		}
		case GZ_AASHIFTX: {
			this->xshiftScalar = *static_cast<float*>(valueList[i]);
			break;
		}
		case GZ_AASHIFTY: {
			this->yshiftScalar = *static_cast<float*>(valueList[i]);
			break;
		}

		default:
			break;
		}

	}
	return GZ_SUCCESS;
}

int GzRender::GzPutTriangle(int numParts, GzToken *nameList, GzPointer *valueList)
/* numParts - how many names and values */
{
/* HW 2.2
-- Pass in a triangle description with tokens and values corresponding to
      GZ_NULL_TOKEN:		do nothing - no values
      GZ_POSITION:		3 vert positions in model space
-- Return error code
*/
/*
-- Xform positions of verts using matrix on top of stack 
-- Clip - just discard any triangle with any vert(s) behind view plane 
		- optional: test for triangles with all three verts off-screen (trivial frustum cull)
-- invoke triangle rasterizer  
*/
	int status = GZ_SUCCESS;
	bool skip = false;
	GzCoord vertListScreen[3];
	GzCoord vertListImg[3];
	GzCoord normalListImg[3];
	GzTextureIndex uvList[3];
	int interpMode = -1;
	int postSortIndex[3]{ 0,1,2 };
	for (unsigned int i = 0; i < numParts; ++i)
	{
		GzToken token{ nameList[i] };
		switch (token)
		{
		case GZ_POSITION:
		{
			GzCoord* vertListModel{ static_cast<GzCoord*>(valueList[i]) };
			RenderUtil::ySort(vertListModel, postSortIndex);
			for (unsigned int j = 0; j < 3; ++j)
			{


				GzCoord4 tri{ vertListModel[j][X],
					vertListModel[j][Y],
				vertListModel[j][Z],
				1.0f };

				// m -> w
				RenderUtil::coordTransform(4, Ximage[matlevel], tri);

				//w -> i
				RenderUtil::coordTransform(4, this->m_camera.Xiw, tri);
				vertListImg[j][X] = tri[X];
				vertListImg[j][Y] = tri[Y];
				vertListImg[j][Z] = tri[Z];
				if (tri[Z] < 0.0f) {
					skip = true;
					break;
				}

				//i -> p
				RenderUtil::coordTransform(4, this->m_camera.Xpi, tri);
				tri[X] = tri[X] / tri[W];
				tri[Y] = tri[Y] / tri[W];
				tri[Z] = tri[Z] / tri[W];
				tri[W] = tri[W] / tri[W];

				//p -> s
				RenderUtil::coordTransform(4, this->Xsp, tri);

				vertListScreen[j][X] = tri[X];
				vertListScreen[j][Y] = tri[Y];
				vertListScreen[j][Z] = tri[Z];

			}

			break;
		}
		case GZ_NORMALS: {
			GzCoord* normalListModel{ static_cast<GzCoord*>(valueList[i]) };
			for (unsigned int j = 0; j < 3; ++j)
			{
				GzCoord4 norm{ normalListModel[postSortIndex[j]][X],
					normalListModel[postSortIndex[j]][Y],
					normalListModel[postSortIndex[j]][Z],
					0.0f
				};

				GzMatrix rotMat;
				RenderUtil::getUrinaryRotation(Ximage[matlevel], rotMat);
				// m -> w
				RenderUtil::coordTransform(4, rotMat, norm);
				//w -> i
				RenderUtil::coordTransform(4, this->m_camera.Xiw, norm);

				normalListImg[j][X] = norm[X];
				normalListImg[j][Y] = norm[Y];
				normalListImg[j][Z] = norm[Z];
				RenderUtil::normalize(normalListImg[j]);

			}

			break;
		}
		case GZ_TEXTURE_INDEX: {
			GzTextureIndex* inUVList{ static_cast<GzTextureIndex*>(valueList[i]) };
			for (unsigned int j = 0; j < 3; ++j)
			{

				uvList[j][U] = inUVList[postSortIndex[j]][U];
				uvList[j][V] = inUVList[postSortIndex[j]][V];
			}
			break;
		}

		default:
			break;
		}

	}
	GzColor vColorList[3];
	GzCoord eye{ 0.0f,0.0f,-1.0f };
	GzColor white{ 1.0f,1.0f,1.0f };
	for (unsigned int vi = 0; vi < 3; ++vi)
	{
		if (this->tex_fun != nullptr) {
			shadeSurface(white, white, white, this->spec, normalListImg[vi], eye, this->lights, this->ambientlight, vColorList[vi]);
		}
		else {
			shadeSurface(this->Ka, this->Ks, this->Kd, this->spec, normalListImg[vi], eye, this->lights, this->ambientlight, vColorList[vi]);
		}
		
	}


	if (!skip) {
		status |= rasterize(vertListScreen, normalListImg,uvList, interpMode, vColorList);
	}


	return status;
}

int GzRender::rasterize(GzCoord* vertList, GzCoord* normalList, GzTextureIndex* uvList, int interpMode, GzColor* vColorList) {
	int status = GZ_SUCCESS;
	float minX, maxX, minY, maxY;
	RenderUtil::getBoundingbox(vertList, minX, maxX, minY, maxY);
	float a, b, c, d;
	status |= RenderUtil::getPlanecoefficient(vertList, a, b, c, d);
	int test{ RenderUtil::LEE(vertList[0],vertList[2],vertList[1][X],vertList[1][Y]) };
	GzCoord* v1;
	GzCoord* v2;
	GzCoord* v3;
	GzCoord norm;

	if (test == L) {
		v1 = &(vertList[2]);
		v2 = &(vertList[1]);
		v3 = &(vertList[0]);
	}
	else {
		v1 = &(vertList[2]);
		v2 = &(vertList[0]);
		v3 = &(vertList[1]);
	}

	float iz;
	GzIntensity red, green, blue, alpha;
	GzDepth z;

	

	//// calculate normal virtual space
	GzCoord normalSpacePointList[3];

	GzCoord xList{ normalList[0][X], normalList[1][X], normalList[2][X] };
	GzCoord yList{ normalList[0][Y], normalList[1][Y], normalList[2][Y] };
	GzCoord zList{ normalList[0][Z], normalList[1][Z], normalList[2][Z] };

	RenderUtil::getVirtualSpacePoint(vertList, xList, normalSpacePointList);
	float nax, nbx, ncx, ndx;
	status |= RenderUtil::getPlanecoefficient(normalSpacePointList, nax, nbx, ncx, ndx);

	RenderUtil::getVirtualSpacePoint(vertList, yList, normalSpacePointList);
	float nay, nby, ncy, ndy;
	status |= RenderUtil::getPlanecoefficient(normalSpacePointList, nay, nby, ncy, ndy);

	RenderUtil::getVirtualSpacePoint(vertList, zList, normalSpacePointList);
	float naz, nbz, ncz, ndz;
	status |= RenderUtil::getPlanecoefficient(normalSpacePointList, naz, nbz, ncz, ndz);

	//// calculate color virtual space
	GzColor colorSpacePointList[3];
	float car, cbr, ccr, cdr;
	float cag, cbg, ccg, cdg;
	float cab, cbb, ccb, cdb;
	if (this->interp_mode == GZ_COLOR) {
		GzColor rList{ vColorList[0][RED],vColorList[1][RED] ,vColorList[2][RED] };
		GzColor gList{ vColorList[0][GREEN],vColorList[1][GREEN] ,vColorList[2][GREEN] };
		GzColor bList{ vColorList[0][BLUE],vColorList[1][BLUE] ,vColorList[2][BLUE] };

		RenderUtil::getVirtualSpacePoint(vertList, rList, colorSpacePointList);

		status |= RenderUtil::getPlanecoefficient(colorSpacePointList, car, cbr, ccr, cdr);

		RenderUtil::getVirtualSpacePoint(vertList, gList, colorSpacePointList);

		status |= RenderUtil::getPlanecoefficient(colorSpacePointList, cag, cbg, ccg, cdg);

		RenderUtil::getVirtualSpacePoint(vertList, bList, colorSpacePointList);

		status |= RenderUtil::getPlanecoefficient(colorSpacePointList, cab, cbb, ccb, cdb);
	}

	////  calculate uv virtual space
	// list of values used for warping uv
	float vzList[3];
	vzList[0] = (vertList[0][Z] / ((float)INT_MAX - vertList[0][Z]))+1.0f;
	vzList[1] = (vertList[1][Z] / ((float)INT_MAX - vertList[1][Z])) + 1.0f;
	vzList[2] =( vertList[2][Z] / ((float)INT_MAX - vertList[2][Z]))+1.0f;
	GzCoord uvSpacePointList[3];
	float au, bu, cu, du;
	float av, bv, cv, dv;
	GzCoord uList{ uvList[0][U]/vzList[0], uvList[1][U] /vzList[1], uvList[2][U] /vzList[2]  };
	GzCoord vList{ uvList[0][V] /vzList[0], uvList[1][V] / vzList[1], uvList[2][V] / vzList[2] };

	//GzCoord uList{ uvList[0][U], uvList[1][U] , uvList[2][U] };
	//GzCoord vList{ uvList[0][V], uvList[1][V], uvList[2][V]};

	RenderUtil::getVirtualSpacePoint(vertList, uList, uvSpacePointList);
	status |= RenderUtil::getPlanecoefficient(uvSpacePointList, au, bu, cu, du);
	RenderUtil::getVirtualSpacePoint(vertList, vList, uvSpacePointList);
	status |= RenderUtil::getPlanecoefficient(uvSpacePointList, av, bv, cv, dv);

	float ir, ig, ib;

	float iNormalX, iNormalY, iNormalZ;

	float iu, iv;

	GzCoord eye{ 0.0f,0.0f,-1.0f };
	for (int x = std::floor(minX); x < std::ceil(maxX)+1; x++)
	{
		for (int y = std::floor(minY); y < std::ceil(maxY)+1; y++)
		{
			float shiftedX = this->xshiftScalar+(float)x;
			float shiftedY = this->yshiftScalar+(float)y;

			if (RenderUtil::boundByTri((*v1), (*v2), (*v3), shiftedX, shiftedY)) {
				// interpolate Z
				RenderUtil::interpolate3d((float)a, (float)b, (float)c, (float)d, shiftedX,shiftedY, iz);
				if (this->interp_mode == GZ_COLOR) {
					// interpolate Color
					RenderUtil::interpolate3d((float)car, (float)cbr, (float)ccr, (float)cdr, shiftedX, shiftedY, ir);
					RenderUtil::interpolate3d((float)cag, (float)cbg, (float)ccg, (float)cdg, shiftedX, shiftedY, ig);
					RenderUtil::interpolate3d((float)cab, (float)cbb, (float)ccb, (float)cdb, shiftedX, shiftedY, ib);
				}

				// interpolate Normal
				RenderUtil::interpolate3d((float)nax, (float)nbx, (float)ncx, (float)ndx, shiftedX, shiftedY, iNormalX);
				RenderUtil::interpolate3d((float)nay, (float)nby, (float)ncy, (float)ndy, shiftedX, shiftedY, iNormalY);
				RenderUtil::interpolate3d((float)naz, (float)nbz, (float)ncz, (float)ndz, shiftedX, shiftedY, iNormalZ);
				GzCoord interpolatedNormal{ iNormalX,iNormalY,iNormalZ };
				RenderUtil::normalize(interpolatedNormal);

				// interpolate UV
				RenderUtil::interpolate3d((float)au, (float)bu, (float)cu, (float)du, shiftedX, shiftedY, iu);
				RenderUtil::interpolate3d((float)av, (float)bv, (float)cv, (float)dv, shiftedX, shiftedY, iv);

				GzGet(shiftedX, shiftedY, &red, &green, &blue, &alpha, &z);
				if (iz < (float)z) {
					GzColor diffuseCol{0.0f,0.0f,0.0f};
					float viz = (iz / ((float)INT_MAX - iz))+1.0f;
					if (this->tex_fun != nullptr) {
						this->tex_fun(iu*viz, iv*viz, diffuseCol);
						//this->tex_fun(iu, iv, diffuseCol);
					}
					if (this->interp_mode == GZ_COLOR) {
						if (this->tex_fun != nullptr) {
							GzPut(shiftedX, shiftedY, ctoi(ir*diffuseCol[RED]), ctoi(ig*diffuseCol[GREEN]), ctoi(ib*diffuseCol[BLUE]), 1, RenderUtil::roundToInt(iz));
						}
						else {
							GzPut(shiftedX, shiftedY, ctoi(ir), ctoi(ig), ctoi(ib), 1, RenderUtil::roundToInt(iz));
						}
						
					}
					else if (this->interp_mode == GZ_NORMALS) {
						
						GzColor phongColor;
						if (this->tex_fun != nullptr) {
							shadeSurface(diffuseCol, this->Ks, diffuseCol, this->spec, interpolatedNormal, eye, this->lights, this->ambientlight, phongColor);
						}
						else {
							shadeSurface(this->Ka, this->Ks, this->Kd, this->spec, interpolatedNormal, eye, this->lights, this->ambientlight, phongColor);
						}
						
						GzPut(shiftedX, shiftedY, ctoi(phongColor[RED]), ctoi(phongColor[GREEN]), ctoi(phongColor[BLUE]), 1, RenderUtil::roundToInt(iz));
					}

				}
			}

		}
	}
	return status;
}