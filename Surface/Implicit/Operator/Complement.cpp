/**
 * Implementation of the complement.
 * @file Complement.cpp
 * @date July 17, 2001
 * @author Ed Bachta
 */

#include "Complement.h"

REGISTER_IMPLICIT(Complement, "UnaryOp:Complement");

#ifndef INTERVAL_EVAL_ONLY
/**
 * Negates the value of the implicit.
 */
double Complement::proc(const gmVector3 & x) 
{ 
  return (m_f ? -m_f->proc(x) : 0.0);
}

/** 
 * Negates the gradient of the implicit.
 */
gmVector3 Complement::grad(const gmVector3 & x) 
{ 
  return (m_f ? -m_f->grad(x) : gmVector3());
}

/** Negates the Hessian of the implicit.
 */
gmMatrix3 Complement::hess(const gmVector3 & x) 
{ 
  return (m_f ? -m_f->hess(x) : gmMatrix3());
}
#endif

Intervald Complement::proc(const Box<double>& x) 
{ 
  return (m_f ? -m_f->proc(x) : Intervald(0.0));
}

Box3d Complement::grad(const Box<double>& x) 
{ 
  if (m_f)
    return -(m_f->grad(x));
  else
    return Box3d(0.0);
}

IMatrix3d Complement::hess(const Box<double>& x) 
{ 
  if (m_f)
    return -(m_f->hess(x));
  else
    return IMatrix3d(0.0);
}

/**
 * Evaluate the derivative of proc with respect to q at a given point.
 * This is the same as the underlying function, with all values negated.
 *
 * @param p    The point at which to evaluate.
 * @param q    The derivatives with respect to q vector elements,
 *             on return.  Must have this->qlen() elements.
 */
void Complement::procq(const gmVector3& p, double * q) 
{
  if (!m_f) 
    return;

  m_f->procq(p,q);
  for(unsigned int i=0; i < m_f->qlen(); i++ )
    q[i] = -q[i];
}

