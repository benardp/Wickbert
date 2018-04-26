#ifndef OPENGL_H
#define OPENGL_H

#include "defs.h"

#if defined(_MSC_VER)  // Identifying MS VC++ compilers
         // Special headers for MS compilers  
	#include <wtypes.h>
	#include <wingdi.h>  

#endif


#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif
//#include <glut.h>

#endif