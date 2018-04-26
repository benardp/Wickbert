/**
 * This file contains the bool (boolean) type definition.  This file is
 * necessary only when compiling with a "C" compiler.
 * @file  bool.h
 * @date  15 June 1994
 * @author Ferdi Scheepers
 * @author Stephen F May
 */

#ifndef BOOL_H
#define BOOL_H

// keep "bool" from interfering with predefined "bool"

#ifdef DEFINE_BOOL

/**
 * Define bool as {false, true} if it isn't defined yet.  This is not used
 * when compiling with C++ compilers.
 */
typedef int bool;
enum { false, true };

#endif

#endif

