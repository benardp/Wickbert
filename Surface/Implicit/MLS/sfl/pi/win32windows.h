/* ### CHECK-ME-BEFORE-RELEASE ### */
/*
 * Title:   win32windows.h
 * Created: Spring 2000
 * Authors: Gilbert Baumann, Markus Noga, Tim Weyrich
 *     $Id: win32windows.h,v 1.2 2003/11/19 10:37:31 weyrich Exp $
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
 * $Log: win32windows.h,v $
 * Revision 1.2  2003/11/19 10:37:31  weyrich
 * added Pointshop3D headers
 *
 *
 */

/*! \file   win32windows.h
    \brief  MS Windows header (somewhat massaged).
    \author Gilbert Baumann, Markus L. Noga, Tim Weyrich
*/

#ifndef __win32windows_h__
#define __win32windows_h__

#if defined(_WIN32) && !defined(_WINDOWS_)
/* MS defines a structure field of this name. That would break
 * platform independence.
 *
 * MS objidl.h is not ANSI compilant - fixing PropVariantInit duplicates linker error.
 */
#  define int64 microsoft_int64

#  ifdef __cplusplus
#   define PropVariantInit __inline PropVariantInit
#   pragma warning( push )
#   pragma warning( disable :  4005 4141 )
#  endif

#  include <windows.h>

#  undef  int64

#  ifdef __cplusplus
#   pragma warning( pop )
#   undef  PropVariantInit
#  endif

#endif

#endif /* __win32windows_h__ */
