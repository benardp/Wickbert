//cleangl.h
//This file is copied from the header of glut.h to make this project cross platform.
//Wen Su
//2004-02-03
#ifndef CLEANGL_H
#define CLEANGL_H


#if defined(_WIN32)
#include <windows.h>
#endif
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif


#endif  // CLEANGL_H
