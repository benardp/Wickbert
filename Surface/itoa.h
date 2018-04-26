/*
 *  itoa.h
 *  GPS
 *
 *  Created by Matei Stroila on 11/12/05.
 *  
 */

//itoa is not a standard C++ function. Here it is an implementation for
//compilers that don't know it -ms

#ifndef ITOA_H
#define ITOA_H
/**

* C++ version char* style "itoa":

*/
char* itoa( int value, char* result, int base );

#endif
