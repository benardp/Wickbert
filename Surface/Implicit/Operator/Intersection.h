/**
 * @file Intersection.h
 * Represents the intersection of two implicit surfaces.
 *
 * @author Steve Zelinka
 */

#ifndef _INTERSECTION_H
#define _INTERSECTION_H

#include "RFunction.h"
#if 0
class Intersection : public RFunction
{
  private:
    void init(Implicit*,Implicit*,double);

  public:
    /// Constructors
    Intersection();
    Intersection(Implicit*, Implicit*);
    Intersection(Implicit*,Implicit*,double);

    MAKE_NAME();
};
#endif

class Intersection : public RFunction
{
public:
	Intersection(Implicit *f, Implicit* g, int cont = 1) : RFunction(cont, 1.0)
	{
		m_f = f;
		m_g = g;
	}
	Intersection() : RFunction(1, 1.0)
	{
		m_f = m_f = NULL;
	}
	
	MAKE_NAME();
};
#endif

