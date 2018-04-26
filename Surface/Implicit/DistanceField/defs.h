#ifndef DEFS_H
#define DEFS_H


#ifdef WIN32
#define WIN32_LEAN_AND_MEAN
//#define NOMINMAX
#include <windows.h>
#endif

#ifndef mkdir
#ifdef WIN32
#include <direct.h>
#define mkdir(a,b) _mkdir(a)
#endif
#endif

#include <assert.h>


#endif
