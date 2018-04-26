/* ### CHECK-ME-BEFORE-RELEASE ### */
/*
 * Title:   win16sound.h
 * Created: Spring 2000
 * Authors: Markus Noga
 *     $Id: win16sound.h,v 1.3 2003/11/19 10:37:31 weyrich Exp $
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
 * $Log: win16sound.h,v $
 * Revision 1.3  2003/11/19 10:37:31  weyrich
 * added Pointshop3D headers
 *
 *
 */

#ifndef __win32sound_h__
#define __win32sound_h__

#include <pisystem.h>
#include <pithread.h>
#include <win32windows.h>

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
  void stop() {
    if(isValid())
      waveOutReset(wout);
  }

  //! Pause playing sound.
  void pause() {
    if(isValid())
      waveOutPause(wout);
  }

  //! Resume from a pause.
  void resume() {
    if(isValid())
      waveOutRestart(wout);
  }

  //! Render next sound buffer and enqueue it.
  /*! Invokes renderBuffer(). Must be public for C
      access. Do not call from your code.
   */
  void renderAndPlayBuffer(WAVEHDR *hdr);

  void renderAndPlayNextBuffer();

  //! Generate a multimedia error.
  void genError(const char *msg,MMRESULT res);

protected:
  PiThread *pth;         //!< Sound computation thread

  HWAVEOUT  wout;	 //!< MMSYSTEM wave output handle

  WAVEHDR  *bufheaders;  //!< MMSYSTEM waveform buffer headers
  uint8    *bufdata;	 //!< Waveform buffer data

  uint32    nextBuffer;
  MMRESULT  timerId;
};

#endif // __win32sound_h__
