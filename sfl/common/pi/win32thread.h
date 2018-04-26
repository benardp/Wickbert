/* ### CHECK-ME-BEFORE-RELEASE ### */
/*
 * Title:   win32thread.h
 * Created: Spring 2000
 * Authors: Markus Noga
 *     $Id: win32thread.h,v 1.2 2003/11/19 10:37:31 weyrich Exp $
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
 * $Log: win32thread.h,v $
 * Revision 1.2  2003/11/19 10:37:31  weyrich
 * added Pointshop3D headers
 *
 *
 */

/*! \file   win32thread.c
    \brief  Win32 thread routines.
    \author Markus L. Noga
*/
#ifndef __win32thread_h__
#define __win32thread_h__

#include <win32windows.h>
#include <pitypes.h>


/*! \brief Thread type
     Internal usage: lower 32 bits: thread handle, higher 32 bits: thread id.
*/
typedef uint64 pithread_t;

/*! \brief Build thread type from handle and id. */
#define MAKE_THREAD(h,id)     ((((uint64)(id))<<32)|(uint64)(h))
/*! \brief Get handle from thread type. */
#define GET_HANDLE(th)	      ((HANDLE)(th))
/*! \brief Get id from thread type. */
#define GET_ID(th)	      ((DWORD)((th)>>32))

/*! \brief Thread attribute type */
typedef uint32 pithread_attr_t;


/*! \brief The mutex type. */
typedef HANDLE  pithread_mutex_t;

/*! \brief The mutex attribute type. */
typedef uint32 pithread_mutexattr_t;


/*! \brief The condition type. */
typedef HANDLE  pithread_cond_t;

/*! \brief The condition attribute type. */
typedef uint32 pithread_condattr_t;


/* Include interface definition
 */
#define CREATE_ATT  static PI_INLINE PI_API
#define GENERAL_ATT static PI_INLINE PI_API

#include <pithread-if.h>


/*! \bug Attributes currently ignored.
*/
CREATE_ATT int pithread_create(pithread_t *__thread,
                              const pithread_attr_t *__attr,
                              void *(*__start_routine) (void *),
                              void *__arg) {
  DWORD  id;
  HANDLE h=CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE) __start_routine,
			__arg, 0, &id);

  if(h)
    *__thread=MAKE_THREAD(h,id);

  return h==NULL;
}

GENERAL_ATT pithread_t pithread_self(void) {
  return MAKE_THREAD(GetCurrentThread(), GetCurrentThreadId());
}

GENERAL_ATT int pithread_equal(pithread_t __thread1, pithread_t __thread2) {
  return GET_ID(__thread1)==GET_ID(__thread2);
}

GENERAL_ATT void pithread_exit(void *__retval) {
  ExitThread((DWORD) __retval);
}

GENERAL_ATT int pithread_join(pithread_t __th, void **__thread_return) {
  HANDLE h=GET_HANDLE(__th);
  DWORD res1=WaitForSingleObject(h, INFINITE);
  
  if(res1==WAIT_OBJECT_0)
    return !GetExitCodeThread(h,(DWORD*) __thread_return);

  return -1;
}

GENERAL_ATT int pithread_cancel(pithread_t __thread) {
  return (int) !TerminateThread(GET_HANDLE(__thread), (DWORD) -1);
}


/*! \bug Defaults to recursive mutex. Attributes are ignored. */
GENERAL_ATT int pithread_mutex_init(pithread_mutex_t *__mutex,
                                    const pithread_mutexattr_t *__mutex_attr) {
  HANDLE h=CreateMutex(NULL,false,NULL);
  if(h) {
    *__mutex=h;
    return 0;
  }
  return -1;
}

GENERAL_ATT int pithread_mutex_destroy(pithread_mutex_t *__mutex) {
  return !CloseHandle(*__mutex);
}

GENERAL_ATT int pithread_mutex_trylock(pithread_mutex_t *__mutex) {
  return !(WaitForSingleObject(*__mutex,0)==WAIT_OBJECT_0);
}

GENERAL_ATT int pithread_mutex_lock(pithread_mutex_t *__mutex) {
  return !(WaitForSingleObject(*__mutex,INFINITE)==WAIT_OBJECT_0);
}

GENERAL_ATT int pithread_mutex_unlock(pithread_mutex_t *__mutex) {
  return !ReleaseMutex(*__mutex);
}






/*! \bug Win32 events don't compare to POSIX conditions. */
GENERAL_ATT int pithread_cond_init(pithread_cond_t *__cond,
                                   const pithread_condattr_t *__cond_attr) {
  HANDLE h=CreateEvent(NULL, 
                       FALSE /* manual reset? */,
		       FALSE /* initially set? */,
		       NULL);
  
  if(h) {
    *__cond=h;
    return 0;
  }
  return -1;
}

/*! \bug Win32 events don't compare to POSIX conditions. */
GENERAL_ATT int pithread_cond_destroy(pithread_cond_t *__cond) {
  return !CloseHandle(*__cond);
}

/*! \bug Win32 events don't compare to POSIX conditions. */
GENERAL_ATT int pithread_cond_signal(pithread_cond_t *__cond) {
  return !SetEvent(*__cond);
}

/*! \bug Unimplemented. */
GENERAL_ATT int pithread_cond_broadcast(pithread_cond_t *__cond) {
  return -1;
}

/*! \bug Unimplemented. */
GENERAL_ATT int pithread_cond_wait(pithread_cond_t *__cond,
                                   pithread_mutex_t *__mutex) {
  return -1;
}


GENERAL_ATT void pithread_yield() {
  Yield();
}

#endif /* __win32thread_h__ */
