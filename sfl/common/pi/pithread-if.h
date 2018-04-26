/* ### CHECK-ME-BEFORE-RELEASE ### */
/*
 * Title:   pithread-if.h
 * Created: Spring 2000
 * Authors: Markus Noga
 *     $Id: pithread-if.h,v 1.2 2003/11/19 10:37:30 weyrich Exp $
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
 * $Log: pithread-if.h,v $
 * Revision 1.2  2003/11/19 10:37:30  weyrich
 * added Pointshop3D headers
 *
 *
 */

/*! \file   pithread-if_h
    \brief  Thread class interface
    \author Markus L. Noga

    Only to be included from the platform implementation headers.
*/

#ifndef __pithread_if_h__
#define __pithread_if_h__

/*! \brief Create a thread.
    \param __thread          storage area for created thread handle.
    \param __attr            thread attributes (or default attributes if NULL)
    \param __start_routine   function to call in thread
    \param __arg             argument to pass to __start_routine
    \return                  0 on success, else error.

    Threads are started with canceltype PTHREAD_CANCEL_ASYNCHRONOUS.
*/
CREATE_ATT int pithread_create(pithread_t *__thread,
                           const pithread_attr_t *__attr,
                           void *(*__start_routine) (void *),
                           void *__arg);

/*! \brief Obtain the identifier of the current thread.
 */
GENERAL_ATT pithread_t pithread_self(void);

/*! \brief Compare two thread identifiers.
    \return nonzero if they refer to the same thread, else zero
 */
GENERAL_ATT int pithread_equal(pithread_t __thread1, pithread_t __thread2);

/*! \brief Terminate calling thread.
    \param __retval  Thread return value
 */
GENERAL_ATT void pithread_exit(void *__retval) PI_NORETURN;

/*! \brief Make calling thread wait for termination of another thread.
    \param __th            thread to wait for.
    \param __thread_return storage area for thread return value or NULL.
    \return                0 on success, else error.
 */
GENERAL_ATT int pithread_join(pithread_t __th, void **__thread_return);

/*! \brief Cancel given thread immediately or at the next possibility.
    \return 0 on success, else error.
 */
GENERAL_ATT int pithread_cancel(pithread_t __thread);



/*! \brief Initialize a mutex.
    \param __mutex      Mutex to initialize
    \param __mutex_attr Attributes to use (default values if NULL).
    \return 0 on success, else error.
 */
GENERAL_ATT int pithread_mutex_init(pithread_mutex_t *__mutex,
                                    const pithread_mutexattr_t *__mutex_attr);

/*! \brief Destroy mutex.
    \param __mutex The mutex.
    \return 0 on success, else error.
 */
GENERAL_ATT int pithread_mutex_destroy(pithread_mutex_t *__mutex);

/*! \brief Try to lock mutex (non-blocking).
    \param __mutex The mutex.
    \return 0 on success, else error.
*/
GENERAL_ATT int pithread_mutex_trylock(pithread_mutex_t *__mutex);

/*! \brief Wait until lock for mutex becomes available and lock it.
    \param __mutex The mutex.
    \return 0 on success, else error.
*/
GENERAL_ATT int pithread_mutex_lock(pithread_mutex_t *__mutex);

/*! \brief Unlock mutex.
    \param __mutex The mutex.
    \return 0 on success, else error.
*/
GENERAL_ATT int pithread_mutex_unlock(pithread_mutex_t *__mutex);



/*! \brief Initialize condition variable.
    \param __cond      condition variable to be initialized
    \param __cond_attr attributes to use (default values if NULL)
    \return 0 on success, else error.
*/
GENERAL_ATT int pithread_cond_init(pithread_cond_t *__cond,
                                   const pithread_condattr_t *__cond_attr);

/*! \brief Destroy condition variable.
    \param __cond  the condition variable
    \return 0 on success, else error.
*/
GENERAL_ATT int pithread_cond_destroy(pithread_cond_t *__cond);

/*! \brief Wake up one thread waiting for condition variable.
    \param __cond  the condition variable
    \return 0 on success, else error.
*/
GENERAL_ATT int pithread_cond_signal(pithread_cond_t *__cond);

/*! \brief Wake up all threads waiting for condition variable.
    \param __cond  the condition variable
    \return 0 on success, else error.

    \bug currently unused.
*/
GENERAL_ATT int pithread_cond_broadcast(pithread_cond_t *__cond);

/*! \brief Wait for condition variable to be signaled or broadcast.
    \param __cond  The condition variable.
    \param __mutex The mutex (assumed to be locked before).
    \return 0 on success, else error.
*/
GENERAL_ATT int pithread_cond_wait(pithread_cond_t *__cond,
                                   pithread_mutex_t *__mutex);


/*! Yield the remainder of the current timeslice to other threads/processes. */
GENERAL_ATT void pithread_yield();

#endif /* __pithread_if_h__ */
