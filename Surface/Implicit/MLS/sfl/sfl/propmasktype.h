// ### CHECK-ME-BEFORE-RELEASE ###
//
// Title:   propmasktype.h
// Created: Sun Nov 16 16:46:05 2003
// Authors: Tim Weyrich <weyrich@inf.ethz.ch>
//
// Copyright (c) 2001--2003, Computer Graphics Lab, ETH Zurich
//
// This file is part of the Pointshop3D system.
// See http://www.pointshop3d.com/ for more information.
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General
// Public License along with this library; if not, write to the
// Free Software Foundation, Inc., 59 Temple Place, Suite 330,
// Boston, MA 02111-1307 USA
//
// Contact info@pointshop3d.com if any conditions of this
// licensing are not clear to you.
//
// ---------------------------------------------------------------
//
// $Log: propmasktype.h,v $
// Revision 1.9  2003/11/21 10:48:02  weyrich
// hopefully nuked g++-3 bug completely...
//
// Revision 1.8  2003/11/20 18:03:43  weyrich
// take this...
//
// Revision 1.7  2003/11/20 15:24:51  weyrich
// fix: argument list termination PropMaskType::setBits()
//
// Revision 1.6  2003/11/20 13:08:04  weyrich
// moved PropMaskType into sfl namespace
//
// Revision 1.5  2003/11/20 12:59:18  weyrich
// fix: set() and clear() resize property bitvector
//
// Revision 1.4  2003/11/19 13:04:02  weyrich
// added Pointshop3D headers
//
// Revision 1.3  2003/11/19 09:50:49  weyrich
// fix: g++ optimization bug
//
// Revision 1.2  2003/11/17 15:07:56  weyrich
// removed g++ warning
//
// Revision 1.1  2003/11/17 15:05:52  weyrich
// new API allowing for more than 32 surfel properties
//

#ifndef __PROPMASKTYPE_H__
#define __PROPMASKTYPE_H__

#include <stdio.h>
#include <string.h>

#include "sfl.h"

namespace sfl
{
  class SFL_API PropMaskType
  {
    int  numBits;

    uint32  *data;

    void  resize(int newNumBits)
    {
      const int  n = (numBits+31) / 32;
      const int  n_new = (newNumBits+31) / 32;

      const int  oldNumBits = numBits;

      if (n_new != n)
        {
          uint32  *newData = new uint32[n_new];
          if (data)
            {
              const int  n_min = n < n_new ? n : n_new;
              memcpy(newData, data, n_min*sizeof(uint32));
              delete[] data;
            }
          data = newData;
        }

      numBits = newNumBits;

      // clear freed or newly allocated bits, respectively
      // (the latter is important to keep comparisions and
      // isZero() fast)
      {
        int  i;
        // clear newly allocated bits:
        for (i=oldNumBits; i<numBits; i++) // ### FIXME: potentially slow
          data[i/32] &= ~(1 << (i & 31));
        // clear remaining bits in last uin32:
        if (0 != (numBits & 31))
          data[n_new-1] &= (0xffffffff >> (32 - (numBits & 31)));
      }
    }

  public:
    PropMaskType()
      : numBits(0),
        data(NULL)
    {
    }

    PropMaskType(const PropMaskType &mask)
      : numBits(mask.numBits)
    {
      const int  n = (numBits+31) / 32;
      data = new uint32[n];
      memcpy(data, mask.data, n*sizeof(uint32));
    }

    PropMaskType(int _numBits, const uint32 *_data)
      : numBits(_numBits)
    {
      const int  n = (numBits+31) / 32;
      data = new uint32[n];
      memcpy(data, _data, n*sizeof(uint32));
    }

    PropMaskType(uint32 oldMask)
    {
      if (oldMask & (1<<31))
        { // ### FIXME: PropMaskType doesn't support "sign bit" yet, so
          // ### let's hack the sign bit extension:
          numBits = 512;
          const int  n = (numBits+31) / 32;
          data = new uint32[n];
          memset(data, 0xff, n*sizeof(uint32));
          *data = oldMask;
          // ### Actually, PropMaskType will never support a sign bit.
          // ### This was only included to keep compatible with code
          // ### that passed '-1' to assign 'all possible properties'.
        }
      else
        {
          numBits = 32;
          data = new uint32[1];
          *data = oldMask;
        }
    }

    ~PropMaskType()
    {
      if (data)
        delete[] data;
    }

    PI_INLINE PropMaskType  &operator=(const PropMaskType &mask)
    {
      const int  nOld = (numBits+31) / 32;
      numBits = mask.numBits;
      const int  n = (numBits+31) / 32;

      if (n != nOld)
        {
          delete[] data;
          data = new uint32[n];
        }

      memcpy(data, mask.data, n*sizeof(uint32));

      return *this;
    }

    const uint32 *getData(void) const { return data; }
    int           getNumBits(void) const { return numBits; }

    PI_INLINE uint32  getUInt32(int idx) const
    {
      const int  n = (numBits+31) / 32;
      if (idx >=0 && idx < n)
        return data[idx];
      else
        return 0;
    }

    PI_INLINE bool  isZero(void) const
    {
      // Relies on unused bits of the last uint32 to be zeroed...
      const int  n = (numBits+31) / 32;
      for (int i=0; i<n; i++)
        if (data[i])
          return false;
      return true;
    }

    PI_INLINE bool  isEqual(const PropMaskType &mask) const
    {
      const int  n = (numBits+31) / 32;
      const int  mask_n = (mask.numBits+31) / 32;

      const int  nMin = n < mask_n ? n : mask_n;

      int  i;

      for (i=0; i<nMin; i++)
        if (data[i] != mask.data[i])
          return false;

      for (; i<n; i++)
        if (data[i])
          return false;

      for (; i<mask_n; i++)
        if (mask.data[i])
          return false;

      return true;
    }

    PI_INLINE bool  operator==(const PropMaskType &mask) const
    {
      return isEqual(mask);
    }

    PI_INLINE bool  operator!=(const PropMaskType &mask) const
    {
      return !isEqual(mask);
    }

    PI_INLINE PropMaskType  &operator&(const PropMaskType &mask)
    {
      const int  n = (numBits+31) / 32;
      const int  mask_n = (mask.numBits+31) / 32;

      const int  nMin = n < mask_n ? n : mask_n;

      int  i;

      for (i=0; i<nMin; i++)
        data[i] &= mask.data[i];

      for (; i<n; i++)
        data[i] = 0;

      return *this;
    }

    PI_INLINE PropMaskType  &operator|(const PropMaskType &mask)
    {
      const int  n = (numBits+31) / 32;
      const int  mask_n = (mask.numBits+31) / 32;

      const int  nMin = n < mask_n ? n : mask_n;

      if (n < mask_n)
        {
          delete[] data;
          data = new uint32[mask_n];

          numBits = mask.numBits;
        }

      int  i;

      for (i=0; i<nMin; i++)
        data[i] |= mask.data[i];

      for (; i<mask_n; i++)
        data[i] = mask.data[i];

      return *this;
    }

    PI_INLINE PropMaskType  &operator^(const PropMaskType &mask)
    {
      const int  n = (numBits+31) / 32;
      const int  mask_n = (mask.numBits+31) / 32;

      const int  nMin = n < mask_n ? n : mask_n;

      if (n < mask_n)
        {
          delete[] data;
          data = new uint32[mask_n];

          numBits = mask.numBits;
        }

      int  i;

      for (i=0; i<nMin; i++)
        data[i] ^= mask.data[i];

      for (; i<mask_n; i++)
        data[i] = mask.data[i];

      return *this;
    }

    PI_INLINE bool  isSet(int idx) const
    {
      if (idx >=0 && idx < numBits)
        return 0 != (data[idx/32] & (1 << (idx & 31)));
      else
        return false;
    }


    PI_INLINE PropMaskType  &set(int idx)
    {
      if (idx >=0)
        {
          if (idx >= numBits)
            resize(idx+1);
          data[idx/32] |= (1 << (idx & 31));
        }
      return *this;
    }

    PI_INLINE PropMaskType  &clear(int idx)
    {
      if (idx >=0)
        {
          if (idx >= numBits)
            resize(idx+1);
          data[idx/32] &= ~(1 << (idx & 31));
        }
      return *this;
    }

    PI_INLINE PropMaskType  &set(int idx, bool value)
    {
      return value ? set(idx) : clear(idx);
    }

    //! Sets a list of bits
    /*! The bit indices are passed as argument list; the last
     *  argument must be '-1'.
     */
    PI_INLINE PropMaskType  &setBits(int idx0, ...)
    {
      va_list        ap;

      va_start(ap, idx0);
      setBitsV(idx0, ap);
      va_end(ap);

      return *this;
    }

    //! Sets a list of bits
    /*! The bit indices are passed as argument list; the last
     *  argument must be '-1'.
     */
    PI_INLINE PropMaskType  &setBitsV(int idx0, va_list ap)
    {
      int  idx;
      for (idx=idx0; idx>=0; idx=va_arg(ap, int))
        set(idx);
      return *this;
    }

    void  dumpHex(FILE *out)
    {
      const int  n = (numBits+31) / 32;
      for (int i=n-1; i>=0; i--)
        fprintf(out, "%08x", (unsigned)data[i]);
    }
  };
};

#endif

// Some Emacs-Hints -- please don't remove:
//
//  Local Variables:
//  mode:C++
//  tab-width:4
//  End:
