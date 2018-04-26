/**
 * This file acts as a wrapper for all of the other libgm header files.
 * When using libgm in your program, you need only include "gm.h".
 * @file  gm.h
 * @date  15 June 1994
 * @author Ferdi Scheepers
 * @author Stephen F May
 */

#ifndef GM_H
#define GM_H

#include "bool.h"
#include "gmConst.h"
#include "gmUtils.h"
#include "gmVec2.h"
#include "gmVec3.h"
#include "gmVec4.h"
#include "gmMat3.h"
#include "gmMat4.h"
#include "gmQuat.h"

/**
 * @mainpage libgm - A Graphics Math Library
 * @section authors Written by
 * Ferdi Scheepers and Stephen F. May <BR>
 * Advanced Computing Center for the Arts and Design <BR>
 * Department of Computer and Information Science <BR>
 * The Ohio State University <BR>
 * Columbus, Ohio, USA 
 * @section description Description
 * libgm provides fundamental data types often required by graphics­related
 * programs.  Included are data types for 2D and 3D vectors, and for 3x3 and
 * 4x4 matrices. In addition, libgm provides a set of constants and utility
 * functions which are useful for graphics programming and a boolean data
 * type bool with boolean constants false and true.  
 * @section constants Constants
 * All constants are double precision, floating point numbers defined to 20
 * digits of precision, although only 6 digits are shown here.
 * <PRE>
 * gm2PI        6.283185...   2*PI 
 * gmDEGTORAD   0.017453...   PI/180 
 * gmE          2.718281...   e 
 * gmEEXPPI    23.140692...   e^PI 
 * gmGOLDEN     1.618033...   golden ratio 
 * gmINVPI      0.318309...   PI^(-1)
 * gmLN10       2.302585...   ln 10
 * gmLN2        0.693147...   ln 2 
 * gmLOG10E     0.434294...   log e 
 * gmLOG2E      1.442695...   lg e 
 * gmPI         3.141592...   PI 
 * gmPIDIV2     1.570796...   PI/2 
 * gmPIDIV4     0.785398...   PI/4 
 * gmRADTODEG  57.295779...   180/PI 
 * gmSQRT2      1.414213...   sqrt(2) 
 * gmSQRT2PI    2.506628...   sqrt(2)*PI
 * gmSQRT3      1.732050...   sqrt(3) 
 * gmSQRT10     3.162277...   sqrt(10)
 * gmSQRTE      1.648721...   sqrt(e)
 * gmSQRTHALF   0.707106...   sqrt(0.5)
 * gmSQRTLN2    0.832554...   sqrt(ln 2)
 * gmSQRTPI     1.772453...   sqrt(PI)
 * gmEPSILON    1.0e­10       next double > 0 
 * gmGOOGOL     1.0e50        large double 
 * </PRE>
 * @section matrices Matrices
 * Matrix elements are stored in row­major order and can be accessed or
 * changed by using two zero­indexed subscript operators.  Note that the
 * adjoint matrix M* = (1/|M|)M^(-1) is used to determine the inverse of a
 * non-singular matrix M^1. 
 * @section Comments Comments
 * Notationally, we use f or fk (k = 1,2,...) to denote variables of type
 * double, and i or ik (k = 1,2,...) to denote variables of type int. We use
 * v or vi (i = 1,2,..) to denote instances of vectors and M or Mi (i =
 * 1,2,..) to denote instances of matrices.  To use the libgm functions,
 * simply include "gm.h" in your code.  This will include all other header
 * files from the libgm source.  If you do not require matrices in your
 * code, then you do not need to add the libgm++ library to your project
 * since all vector methods are written using "inline" and are completely
 * contained in the header files.  
 */

#endif 

