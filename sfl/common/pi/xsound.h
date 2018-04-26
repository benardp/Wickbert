/* ### CHECK-ME-BEFORE-RELEASE ### */
/*
 * Title:   xsound.h
 * Created: Spring 2000
 * Authors: Markus Noga
 *     $Id: xsound.h,v 1.2 2003/11/19 10:37:31 weyrich Exp $
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
 * $Log: xsound.h,v $
 * Revision 1.2  2003/11/19 10:37:31  weyrich
 * added Pointshop3D headers
 *
 *
 */

#ifndef __xsound_h__
#define __xsound_h__

#include <pithread.h>
#include <pisystem.h>

//! Bind platform-independent interface to X implementation.
#define XSound  PiSound

//! X thread class implementation.
class XSound : public IfSound {
public:
  //! Constructor.
  /*! \param device   The device file name, (void*)0 -> default.
      \param rate     Sampling rate in Hertz.
      \param bits     Width of a channel in bits, 8 or 16.
      \param channels Number of channels (1=mono, 2=stereo).
      \param bufSize  Size of a playback buffer in samples.
      \param buffers  Number of playback buffers to allocate (usually 2).
   */
  XSound(const char *device,
         uint32 rate, uint32 bits, uint32 channels,
         uint32 bufsize, uint32 buffers);

  ~XSound();

  //! Start playing sound.
  void play();

  //! Stop playing sound.
  void stop();

  //! Render next sound buffer and enqueue it.
  /*! Invokes renderBuffer(). Must be public for C
      access. Do not call from your code.
   */
  void renderAndPlayBuffer();

  //! Buffer-flipping sync.
  /*! Must be public for C access. Do not call from your code.
   */
  void syncBuffer();

protected:
  //! set select(...) I/O trigger
  void setEnable(bool enable);

  //! Return number of bytes that can be written w/o blocking.
  uint32 getNonBlocking();

  //! Reset sound system.
  void reset();

  //! Generate an error message from msg and C error number.
  void genError(const char *msg,int nr);

protected:
  pithread_t pth;   //!< sound I/O thread.

  int fd;           //!< sound file descriptor

  uint8 *bufdata;   //!< sound data buffer
};

#endif // __xsound_h__
