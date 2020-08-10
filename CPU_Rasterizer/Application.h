// Application.h: interface for the Application class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_APPLICATION_H__3387B79A_B69F_491D_B782_81D9CAFAAB0F__INCLUDED_)
#define AFX_APPLICATION_H__3387B79A_B69F_491D_B782_81D9CAFAAB0F__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Gz.h"
#include "rend.h"
#include <vector>

class Application  
{
public:
	Application();
	virtual ~Application();
	
public:
	GzRender*  m_pRender;		// the renderer
	GzInput*   m_pUserInput;
	char* m_pFrameBuffer;	// Frame Buffer
	int   m_nWidth;			// width of Frame Buffer
	int   m_nHeight;		// height of Frame Buffer
	std::vector<GzRender*> rtList;		// list of intermediate render target for anti-aliazing
	const float AAFilter[AAKERNEL_SIZE][3]{
	{-0.52f,0.38f,0.128f},{0.41f,0.56f,0.119f},{0.27f,0.08f,0.294f},
	{-0.17f,-0.29f,0.249},{0.58f,-0.55f,0.104},{-0.31f,-0.71f,0.106f}
	}; // kernal for anti-aliazing filter

	virtual int Render()=0; // Pass user input data and call renderer
};

#endif // !defined(AFX_APPLICATION_H__3387B79A_B69F_491D_B782_81D9CAFAAB0F__INCLUDED_)
