/*
@file nv_pbuffer.h
extracted from simple shadow map example from nvidia
for more inforation see PixelBuffers.pdf
*/
#ifndef NV_PBUFFER
#define NV_PBUFFER

#include <gl/glew.h>
#ifdef _WIN32
#include <gl/wglew.h>
#endif

// Define a PBuffer object.
class nv_pbuffer
{

public:

	HPBUFFERARB  hpbuffer;      // Handle to a pbuffer.
    HDC          hdc;           // Handle to a device context.
    HGLRC        hglrc;         // Handle to a GL rendering context.
    int          width;         // Width of the pbuffer
    int          height;        // Height of the pbuffer

	bool		 onecontext;	// A flag to tell us if the pbuffer is sharing 
								// an OpenGL context with our rendering context.

	// constructor and destructor
	nv_pbuffer(int width, int height, int pass);
	~nv_pbuffer();

	// methods
	void testExtensions();
	void wglGetLastError();

};

#endif