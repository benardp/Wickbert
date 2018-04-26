/* ### CHECK-ME-BEFORE-RELEASE ### */
/*
 * Title:   pisystem.h
 * Created: Spring 2000
 * Authors: Gilbert Baumann, Markus Noga, Tim Weyrich
 *     $Id: pisystem.h,v 1.8 2003/11/19 10:37:30 weyrich Exp $
 *
 * Copyright (c) 2000--2003
 * Gilbert Baumann <unk6@stud.uni-karlsruhe.de>
 * Markus Noga <markus@noga.de>
 * Tim Weyrich <weyrich@inf.ethz.ch>
 *
 * This file is part of the Pointshop3D system.
 * See http: *www.pointshop3d.com/ for more information.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330,
 * Boston, MA 02111-1307 USA
 *
 * Contact info@pointshop3d.com if any conditions of this
 * licensing are not clear to you.
 *
 * ---------------------------------------------------------------
 *
 * $Log: pisystem.h,v $
 * Revision 1.8  2003/11/19 10:37:30  weyrich
 * added Pointshop3D headers
 *
 *
 */

/*! \file   pisystem_h
    \brief  Platform-independent system routines.
    \author Gilbert Baumann, Markus L. Noga, Tim Weyrich
*/

#if defined( _WIN32 ) && !defined( __pisocket_h__ )
#  include "pisocket.h"  /* to prevent winsock/winsock2 botch */
#endif

#ifndef __pisystem_h__
#define __pisystem_h__

#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>

#ifdef _WIN32
#  include <malloc.h>
#  define  alloca  _alloca
#  include <io.h>
#  include <process.h>
#else
#  include <unistd.h>
#  include <alloca.h>
#endif

#include "pitypes.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _WIN32
/* By default, Win32 only defines versions of these functions with an
   underscore prefix. We alleviate that.
 */

#  define O_BINARY  _O_BINARY
#  define O_CREAT   _O_CREAT
#  define O_WRONLY  _O_WRONLY
#  define O_TRUNC   _O_TRUNC
#  define O_EXCL    _O_EXCL
typedef int mode_t;

typedef long off_t;

#  define  fileno  _fileno

_CRTIMP int __cdecl access(const char *, int);
_CRTIMP int __cdecl chmod(const char *, int);
_CRTIMP int __cdecl chsize(int, long);
_CRTIMP int __cdecl close(int);
_CRTIMP int __cdecl creat(const char *, int);
_CRTIMP int __cdecl dup(int);
_CRTIMP int __cdecl dup2(int, int);
_CRTIMP int __cdecl eof(int);
_CRTIMP long __cdecl filelength(int);
_CRTIMP int __cdecl isatty(int);
_CRTIMP int __cdecl locking(int, int, long);
_CRTIMP off_t __cdecl lseek(int, off_t, int);
_CRTIMP char * __cdecl mktemp(char *);
_CRTIMP int __cdecl open(const char *, int, ...);
_CRTIMP int __cdecl read(int, void *, unsigned int);
_CRTIMP int __cdecl setmode(int, int);
_CRTIMP int __cdecl sopen(const char *, int, int, ...);
_CRTIMP long __cdecl tell(int);
_CRTIMP int __cdecl umask(int);
_CRTIMP int __cdecl unlink(const char *);
_CRTIMP int __cdecl write(int, const void *, unsigned int);

static PI_INLINE pipe(int *pipes) {
  return _pipe(pipes,512,O_BINARY);
}

#  define popen _popen

#  define fdopen _fdopen

#  define O_RDONLY _O_RDONLY
#  define O_RDWR _O_RDWR

#else
/* Windows needs a rather stupid flag to read in binary. It is undefined under Unix.
 */
#  define  O_BINARY  0

/* Create a new stream connected to a pipe running the given command.  */
  extern FILE *popen (const char *__command, const char *__modes) __PI_THROW;

/* Close a stream opened by popen and return the status of its child.  */
  extern int pclose (FILE *__stream) __PI_THROW;

/* Create a new stream that refers to an existing system file descriptor.  */
  extern FILE *fdopen (int __fd, const char *__modes) __PI_THROW;

#endif /* _WIN32 */

#if defined(_WIN32) || defined (__STRICT_ANSI__)
/* Windows and ANSI-C don't know srandom and random, only srand and rand.
 */
/* Now I found a srandom/random-definition on SGIs, that collides
 * with the following defines. For this reason I just /forbid/ the
 * use of random/srandom. Use rand/srand instead... [tim]
 */
/*
#  define srandom(a)  srand(a)
#  define random()    rand()
*/
#endif

/* FIXME: SGI is defined by my Makefile structur. I should find a
   better way to determine the current platform */
#ifdef SGI
/* Get the pathname of the current working directory,
   and put it in SIZE bytes of BUF.  Returns NULL if the
   directory couldn't be determined or SIZE was too small.
   If successful, returns BUF.  In GNU, if BUF is NULL,
   an array is allocated with `malloc'; the array is SIZE
   bytes long, unless SIZE == 0, in which case it is as
   big as necessary.  */
extern char *pi_getcwd(char *__buf, size_t __size) __PI_THROW;
#else
#  define pi_getcwd getcwd
#endif

#ifdef _WIN32
#  define  PI_PATHSEP  '\\'
#  define  PI_PATHSEP_STR  "\\"
#else
#  define  PI_PATHSEP  '/'
#  define  PI_PATHSEP_STR  "/"
#endif
  
#ifdef __cplusplus
}
#endif


#endif /* __pisystem_h__ */
