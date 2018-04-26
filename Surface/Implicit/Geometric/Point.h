/**
 * Declaration of the geometric point.
 * @file Point.h
 * @date Fall 2000
 * @author Ed Bachta
 * @author Terry Fleury
 */

#ifndef __ASL_POINT_H__
#define __ASL_POINT_H__

#include "Sphere.h"

class Point : public Sphere
{
  public:
    /// Constructors
    Point();
    Point(const gmVector3 & x);

    virtual void _setq(double*);       ///< Sets parameters
    virtual void getq(double* q);      ///< Gets parameters
    virtual unsigned int qlen() { return 3; }   ///< Returns # of parameters
    virtual void getqname(char **qn);  ///< Get parameter names

    MAKE_NAME();
};

#endif

