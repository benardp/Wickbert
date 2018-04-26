//only included if WB_USE_CG is 1
#ifdef WB_USE_CG

// make sure glh_ext is only compiled once
#include "nv_pbuffer.h"
#include <iostream>
using namespace std;

#define REQUIRED_EXTENSIONS "WGL_ARB_extensions_string " \
							"WGL_ARB_pbuffer " \
                            "WGL_NV_render_texture_rectangle " \
                            "WGL_ARB_pixel_format " \
                            "WGL_NV_render_depth_texture " \
							"GL_ARB_multitexture "
						

#define MAX_ATTRIBS     256
#define MAX_PFORMATS    256

/* This routine doesn't seem to find the WGL extensions so it will fail
 * but not bomb out. */
void nv_pbuffer::testExtensions()
{
	cerr << glGetString(GL_EXTENSIONS) << endl;

	if(!glewIsSupported(REQUIRED_EXTENSIONS))
	{
		char *ext = (char *)glGetString(GL_EXTENSIONS);

		cerr << "Some necessary extensions were not supported:" << endl
			 << REQUIRED_EXTENSIONS;
		cerr << "Supported Extensions: " << glGetString(GL_EXTENSIONS) << endl << endl;
		cerr << "Supported Renderer: " << glGetString(GL_RENDERER) << endl << endl;
		//exit(0);
	}
}

nv_pbuffer::~nv_pbuffer()
{

    if ( hpbuffer )
    {
		// Only delete the context if it belongs to us.
        if (onecontext == FALSE)
			wglDeleteContext( hglrc );

        wglReleasePbufferDCARB(hpbuffer, hdc);
		wglGetLastError();


        wglDestroyPbufferARB(hpbuffer);
		wglGetLastError();
        hpbuffer = 0;
   }
}

nv_pbuffer::nv_pbuffer(int width, int height, int pass)
{
    memset(this,0,sizeof(nv_pbuffer));

	HDC hdc = wglGetCurrentDC();
	int pixelformat = GetPixelFormat(hdc);
	HGLRC hglrc = wglGetCurrentContext();
    wglGetLastError();

    // Query for a suitable pixel format based on the specified mode.
    int     format;
    int     pformat[MAX_PFORMATS];
    unsigned int nformats;
	int     iattributes[2*MAX_ATTRIBS];
    float   fattributes[2*MAX_ATTRIBS];
    int     nfattribs = 0;
    int     niattribs = 0;

    // Attribute arrays must be "0" terminated - for simplicity, first
    // just zero-out the array entire, then fill from left to right.
    memset(iattributes,0,sizeof(int)*2*MAX_ATTRIBS);
    memset(fattributes,0,sizeof(float)*2*MAX_ATTRIBS);

    // Since we are trying to create a pbuffer, the pixel format we
    // request (and subsequently use) must be "p-buffer capable".
    iattributes[niattribs  ] = WGL_DRAW_TO_PBUFFER_ARB;
    iattributes[++niattribs] = GL_TRUE;

	// Since we're trying to create a depth texture of the entire
	// screen (which may have an arbitrary aspect ratio) we need
	// a texture rectangle.
	
	if (pass == 1) {
		iattributes[++niattribs] = WGL_BIND_TO_TEXTURE_RECTANGLE_DEPTH_NV;
		iattributes[++niattribs] = GL_TRUE;
	}
	else if (pass == 2) {

		iattributes[++niattribs] = WGL_RED_BITS_ARB;
		iattributes[++niattribs] = 16;

		iattributes[++niattribs] = WGL_GREEN_BITS_ARB;
		iattributes[++niattribs] = 16;

		iattributes[++niattribs] = WGL_BLUE_BITS_ARB;
		iattributes[++niattribs] = 16;

		iattributes[++niattribs] = WGL_ALPHA_BITS_ARB;
		iattributes[++niattribs] = 16;
		
		iattributes[++niattribs] = WGL_FLOAT_COMPONENTS_NV;
		iattributes[++niattribs] = GL_TRUE;

		iattributes[++niattribs] = WGL_BIND_TO_TEXTURE_RECTANGLE_FLOAT_RGBA_NV;
		iattributes[++niattribs] = GL_TRUE;
	}

	PFNWGLCHOOSEPIXELFORMATARBPROC wglChoosePixelFormatARB = 
		(PFNWGLCHOOSEPIXELFORMATARBPROC)wglGetProcAddress("wglChoosePixelFormatARB");

    if ( !wglChoosePixelFormatARB( hdc, 
                                   iattributes, 
                                   fattributes, 
                                   MAX_PFORMATS, 
                                   pformat, 
                                   &nformats ) )
    {
        cerr << "pbuffer creation error:  wglChoosePixelFormatARB() failed.\n";
        exit( -1 );
    }
    wglGetLastError();
	if ( nformats <= 0 )
    {
        cerr << "pbuffer creation error:  Couldn't find a suitable pixel format.\n";
        exit( -1 );
    }


	// Ideally, we'd like to share the same context we're currently using for our
	// rendering window with the pbuffer.  In order to do this, the pixel formats
	// must be the same for both the pbuffer and the rendering window.  So, we check
	// to see if whatever pixel format we're using for our rendering window is in the
	// list of appropriate formats we just created.  If it, we use the existing context
	// and if not, we'll create a new context.

	format = pformat[0];
	
	if (pass == 1) {
		for (unsigned int i = 0; i < MAX_PFORMATS; i++) {
			if (pformat[i] == pixelformat) {
				format = pixelformat;
				onecontext = true;
				break;
			}
		}
	}
	
    // set the pbuffer attributes...
    memset(iattributes,0,sizeof(int)*2*MAX_ATTRIBS);
    niattribs = 0;

    // we need to render to a depth texture
	if (pass == 1) {
		iattributes[niattribs] = WGL_DEPTH_TEXTURE_FORMAT_NV;
		iattributes[++niattribs] = WGL_TEXTURE_DEPTH_COMPONENT_NV;
	}
	else if (pass == 2) {
		iattributes[niattribs] = WGL_TEXTURE_FORMAT_ARB;
		iattributes[++niattribs] = WGL_TEXTURE_FLOAT_RGBA_NV;
	}

	// and we need to render to a rectangle
	iattributes[++niattribs] = WGL_TEXTURE_TARGET_ARB;
    iattributes[++niattribs] = WGL_TEXTURE_RECTANGLE_NV;

	PFNWGLCREATEPBUFFERARBPROC wglCreatePbufferARB = 
		(PFNWGLCREATEPBUFFERARBPROC)wglGetProcAddress("wglCreatePbufferARB");

    // Create the p-buffer.
    this->hpbuffer = wglCreatePbufferARB( hdc, format, width, height, iattributes );
    if ( this->hpbuffer == 0)
    {
        cerr << "pbuffer creation error:  wglCreatePbufferARB() failed\n";
        wglGetLastError();
        return;
    }
    wglGetLastError();

	PFNWGLGETPBUFFERDCARBPROC wglGetPbufferDCARB = 
		(PFNWGLGETPBUFFERDCARBPROC)wglGetProcAddress("wglGetPbufferDCARB");

    // Get the device context.
    this->hdc = wglGetPbufferDCARB( this->hpbuffer );
    if ( this->hdc == 0)
    {
        cerr << "pbuffer creation error:  wglGetPbufferDCARB() failed\n";
        wglGetLastError();
		return;
    }
    wglGetLastError();

	// Create a context for the p-buffer IF we're not using the same context as the
	// current rendering window.
	if (onecontext == true)
		this->hglrc = hglrc;
	else 
		this->hglrc = wglCreateContext( this->hdc );
    
	if ( this->hglrc == 0)
    {
		cerr << "pbuffer creation error:  wglCreateContext() failed\n";
		wglGetLastError();
		return;
    }
    wglGetLastError();

	PFNWGLQUERYPBUFFERARBPROC wglQueryPbufferARB =
		(PFNWGLQUERYPBUFFERARBPROC)wglGetProcAddress("wglQueryPbufferARB");

    // Determine the actual width and height we were able to create.
    wglQueryPbufferARB( this->hpbuffer, WGL_PBUFFER_WIDTH_ARB, &this->width );
    wglQueryPbufferARB( this->hpbuffer, WGL_PBUFFER_HEIGHT_ARB, &this->height );
}

// If we get an error, crash the program!
void nv_pbuffer::wglGetLastError()
{
    DWORD err = GetLastError();
    switch(err)
    {
    case ERROR_INVALID_PIXEL_FORMAT:
        cerr << "Win32 Error:  ERROR_INVALID_PIXEL_FORMAT\n";
		exit(-1);
        break;
    case ERROR_NO_SYSTEM_RESOURCES:
        cerr << "Win32 Error:  ERROR_NO_SYSTEM_RESOURCES\n";
        exit(-1);
		break;
    case ERROR_INVALID_DATA:
        cerr << "Win32 Error:  ERROR_INVALID_DATA\n";
        exit(-1);
		break;
    case ERROR_INVALID_WINDOW_HANDLE:
        cerr << "Win32 Error:  ERROR_INVALID_WINDOW_HANDLE\n";
        exit(-1);
		break;
    case ERROR_RESOURCE_TYPE_NOT_FOUND:
        cerr << "Win32 Error:  ERROR_RESOURCE_TYPE_NOT_FOUND\n";
        exit(-1);
		break;
    case ERROR_SUCCESS:
        // no error
        break;
    default:
        cerr << "Win32 Error:  " << err << endl;
        exit(-1);
		break;
    }
    SetLastError(0);
}
#endif

