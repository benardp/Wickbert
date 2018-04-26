/* ### CHECK-ME-BEFORE-RELEASE ### */
/*
 * Title:   xthread.h
 * Created: Spring 2000
 * Authors: Markus Noga
 *     $Id: xthread.h,v 1.2 2003/11/19 10:37:31 weyrich Exp $
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
 * $Log: xthread.h,v $
 * Revision 1.2  2003/11/19 10:37:31  weyrich
 * added Pointshop3D headers
 *
 *
 */

/*! \file   xthread_h
    \brief  UNIX/X11 implementation of thread functions.
    \author Markus L. Noga

    This file is only to be included from pithread.h on UNIX/X11
*/

#ifndef __xthread_h__
#define __xthread_h__

#include <pthread.h>
#include <pisystem.h>

/*! \brief The thread type. */
typedef pthread_t        pithread_t;

/*! \brief The thread attribute type. */
typedef pthread_attr_t   pithread_attr_t;


/*! \brief The mutex type. */
typedef pthread_mutex_t  pithread_mutex_t;

/*! \brief The mutex attribute type. */
typedef pthread_mutexattr_t pithread_mutexattr_t;


/*! \brief The condition type. */
typedef pthread_cond_t  pithread_cond_t;

/*! \brief The condition attribute type. */
typedef pthread_condattr_t pithread_condattr_t;


/* Include interface definition
 */
#define CREATE_ATT  extern PI_API
#define GENERAL_ATT static PI_INLINE PI_API

#include <pithread-if.h>


GENERAL_ATT pithread_t pithread_self(void) {
  return pthread_self();
}

GENERAL_ATT int pithread_equal(pithread_t __thread1, pithread_t __thread2) {
  return pthread_equal(__thread1,__thread2);
}

GENERAL_ATT void pithread_exit(void *__retval) {
  pthread_exit(__retval);
}

GENERAL_ATT int pithread_join(pithread_t __th, void **__thread_return) {
  return pthread_join(__th,__thread_return);
}

GENERAL_ATT int pithread_cancel(pithread_t __thread) {
  return pthread_cancel(__thread);
}


GENERAL_ATT int pithread_mutex_init(pithread_mutex_t *__mutex,
                                    const pithread_mutexattr_t *__mutex_attr) {
  return pthread_mutex_init(__mutex,__mutex_attr);
}

GENERAL_ATT int pithread_mutex_destroy(pithread_mutex_t *__mutex) {
  return pthread_mutex_destroy(__mutex);
}

GENERAL_ATT int pithread_mutex_trylock(pithread_mutex_t *__mutex) {
  return pthread_mutex_trylock(__mutex);
}

GENERAL_ATT int pithread_mutex_lock(pithread_mutex_t *__mutex) {
  return pthread_mutex_lock(__mutex);
}

GENERAL_ATT int pithread_mutex_unlock(pithread_mutex_t *__mutex) {
  return pthread_mutex_unlock(__mutex);
}



GENERAL_ATT int pithread_cond_init(pithread_cond_t *__cond,
                                   const pithread_condattr_t *__cond_attr) {
  return pthread_cond_init(__cond, __cond_attr);
}

GENERAL_ATT int pithread_cond_destroy(pithread_cond_t *__cond) {
  return pthread_cond_destroy(__cond);
}

GENERAL_ATT int pithread_cond_signal(pithread_cond_t *__cond) {
  return pthread_cond_signal(__cond);
}

GENERAL_ATT int pithread_cond_broadcast(pithread_cond_t *__cond) {
  return pthread_cond_broadcast(__cond);
}

GENERAL_ATT int pithread_cond_wait(pithread_cond_t *__cond,
                                   pithread_mutex_t *__mutex) {
  return pthread_cond_wait(__cond, __mutex);
}

GENERAL_ATT void pithread_yield() {
#ifdef LINUX
  sched_yield();
#else
  sleep(0);
#endif
}

#endif /* __xthread_h__ */
