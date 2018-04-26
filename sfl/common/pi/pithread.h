/* ### CHECK-ME-BEFORE-RELEASE ### */
/*
 * Title:   pithread.h
 * Created: Spring 2000
 * Authors: Markus Noga
 *     $Id: pithread.h,v 1.2 2003/11/19 10:37:30 weyrich Exp $
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
 * $Log: pithread.h,v $
 * Revision 1.2  2003/11/19 10:37:30  weyrich
 * added Pointshop3D headers
 *
 *
 */

/*! \file   pithread_h
    \brief  Thread functions.
    \author Markus L. Noga
*/

#ifndef __pithread_h__
#define __pithread_h__

#ifdef __cplusplus
extern "C" {
#endif

#include <pitypes.h>

#ifdef _WIN32
# include <win32thread.h>
#else
# include <xthread.h>
#endif

#ifdef __cplusplus
}
#endif

#endif /* __pithread_h__ */
