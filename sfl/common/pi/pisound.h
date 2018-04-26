/* ### CHECK-ME-BEFORE-RELEASE ### */
/*
 * Title:   pisound.h
 * Created: Spring 2000
 * Authors: Markus Noga
 *     $Id: pisound.h,v 1.2 2003/11/19 10:37:30 weyrich Exp $
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
 * $Log: pisound.h,v $
 * Revision 1.2  2003/11/19 10:37:30  weyrich
 * added Pointshop3D headers
 *
 *
 */

/*! \file   pisound_h
    \brief  Sound output class.
    \author Markus L. Noga
*/

#ifndef __pisound_h__
#define __pisound_h__

#include "pitypes.h"
#include "pistring.h"
#include "pisystem.h"

// Forward declaration of the platform-independent class.
//
class PiSound;

//! sound class interface.
/*! This is an abstract interface, never to be instantiated.

    Do not use If_____ classes directly. Don't even write
    If_____ in your programs. Instead, use this interface on
    the appropriate Pi_____ class. This applies for class
    methods, too.

    As platforms never change at run time, virtual functions
    aren't justified because of their overhead.
 */
class IfSound {
public:
  //! Constructor.
  /*! \param device   The device file name
      \param rate     Sampling rate in Hertz.
      \param bits     Width of a channel in bits, 8 or 16.
      \param channels Number of channels (1=mono, 2=stereo).
      \param bufSize  Size of a playback buffer in samples
      \param buffers  Number of playback buffers to allocate (usually 2).
   */
  IfSound(const char *device,
          uint32 rate, uint32 bits, uint32 channels,
          uint32 bufsize, uint32 buffers);

  //! Destructor.
  virtual ~IfSound();

  //! Start playing sound.
  void play();

  //! Stop playing sound.
  void stop();

  // Member variable accessors.
  //
  bool        isValid() const       { return valid; }
  const char *getErrMsg() const { return errMsg; }

  const char *getDeviceName() const { return deviceName; }
  uint32      getRate() const       { return rate; }
  uint32      getBits() const       { return bits; }
  uint32      getChannels() const   { return channels; }
  uint32      getSampleSize() const { return sampleSize; }
  uint32      getBufSize() const    { return bufSize; }
  uint32      getBuffers() const    { return buffers; }

  uint32      getBufSizeBytes() const { return bufSize*sampleSize; }

protected:
  //! Render next sound buffer
  /*! To be provided by a subclass. The buffer holds bufSize samples
      of getSampleSize() bytes each.
   */
  virtual void renderBuffer(void *buffer, uint32 bufSize) = 0;

  //! Generate an error message
  void genError(const char *msg);

protected:
  bool   valid;           //!< flag: valid so far.
  char  *errMsg;          //!< sound error message, if any.

  char  *deviceName;      //!< sound device name

  uint32 rate;
  uint32 bits;
  uint32 channels;
  uint32 sampleSize;      //!< Size of a sample in bytes
  uint32 bufSize;
  uint32 buffers;
};


// Bind to implementation for the appropriate platform.
//
#ifdef _WIN32
# include "win32sound.h"
#else
# include "xsound.h"
#endif

#endif // __pisound_h__
