/* ### CHECK-ME-BEFORE-RELEASE ### */
/*
 * Title:   pisocket.h
 * Created: Spring 2000
 * Authors: Markus Noga
 *     $Id: pisocket.h,v 1.13 2003/11/19 10:37:30 weyrich Exp $
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
 * $Log: pisocket.h,v $
 * Revision 1.13  2003/11/19 10:37:30  weyrich
 * added Pointshop3D headers
 *
 *
 */

/*! \file   pisocket_h
    \brief  platform-independent socket services
    \author Markus L. Noga
 */

#ifndef __pisocket_h__
#define __pisocket_h__

#include "pitypes.h"
#include "pisystem.h"

#ifdef _WIN32
/*# include <win32windows.h>*/
# include "winsock.h" //ACC - originally included winsock2.h, caused conflicts
# include "sys/stat.h"

#define  S_IRUSR   _S_IREAD
#define  S_IREAD   _S_IREAD
#define  S_IWUSR   _S_IWRITE
#define  S_IWRITE  _S_IWRITE
#define  S_IXUSR   0100
#define  S_IEXEC   0100
#define  S_IRWXU   (S_IRUSR | S_IWUSR | S_IXUSR)
#define  S_IRGRP   040
#define  S_IWGRP   020

#else
# include <arpa/inet.h> /* Marco Nef <marco@shima.ch> */
# include <netdb.h>
# include <netinet/in.h>
# include <sys/stat.h>
# include <sys/socket.h>
# include <sys/wait.h>
#endif

#if defined( _SGIAPI ) || defined( _WIN32 )
#  define  socklen_t int
#endif

#endif /* __pisocket_h__ */
