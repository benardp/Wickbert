/* ### CHECK-ME-BEFORE-RELEASE ### */
/*
 * Title:   win32sound.h
 * Created: Spring 2000
 * Authors: Markus Noga
 *     $Id: win32sound.h,v 1.3 2003/11/19 10:37:31 weyrich Exp $
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
 * $Log: win32sound.h,v $
 * Revision 1.3  2003/11/19 10:37:31  weyrich
 * added Pointshop3D headers
 *
 *
 */

#ifndef __win32sound_h__
#define __win32sound_h__

#include <pisystem.h>
#include <win32windows.h>
#include <dsound.h>

//! Bind platform-independent interface to Win32 implementation.
#define Win32Sound  PiSound

//! X thread class implementation.
class Win32Sound : public IfSound {
public:
  Win32Sound(const char *device,
             uint32 rate, uint32 bits, uint32 channels,
             uint32 bufsize, uint32 buffers);

  ~Win32Sound();

  //! Start playing sound.
  void play();

  //! Stop playing sound.
  void stop();

  //! Generate a multimedia error.
  void genError(const char *msg,MMRESULT res);

  //! Timer event callback.
  void timerCallback();

protected:
  //! Return write position in samples (not bytes!)
  uint32 getWritePos() {
    uint32 tmp;
    buffer->GetCurrentPosition(NULL,&tmp);
    return tmp/getSampleSize();
  }

  uint32 getPlayPos() {
    uint32 tmp;
    buffer->GetCurrentPosition(&tmp,NULL);
    return tmp/getSampleSize();
  }


  LPDIRECTSOUND	      ds;      //!< Direct sound device

  LPDIRECTSOUNDBUFFER buffer;  //!< Direct sound buffer
  bool                primary; //!< Flag: using primary buffer
  uint32              bufPos;  //!< Playback position in buffer.

  MMRESULT	      timerId; //!< Playback timer
  uint32              timerInt;//!< Timer interval in ms.
};

#endif // __win32sound_h__
