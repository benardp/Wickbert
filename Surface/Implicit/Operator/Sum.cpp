/**
 * @file Sum.cpp
 * Implementation of sum operation
 */
#include "Sum.h"

REGISTER_IMPLICIT(Sum,"BinaryOp:Sum");

#ifndef INTERVAL_EVAL_ONLY
double Sum::proc(const gmVector3 & x) 
{
  return ((m_f ? m_f->proc(x) : 0.0) + 
          (m_g ? m_g->proc(x) : 0.0));
}

gmVector3 Sum::grad(const gmVector3 & x) 
{
  return ((m_f ? m_f->grad(x) : gmVector3()) + 
          (m_g ? m_g->grad(x) : gmVector3())); 
}

gmMatrix3 Sum::hess(const gmVector3 & x) 
{
  return ((m_f ? m_f->hess(x) : gmMatrix3()) + 
          (m_g ? m_g->hess(x) : gmMatrix3())); 
}
#endif

Intervald Sum::proc(const Box<double>&  b) 
{ 
  return ((m_f ? m_f->proc(b) : Intervald(0.0)) + 
          (m_g ? m_g->proc(b) : Intervald(0.0))); 
}

Box3d Sum::grad(const Box<double>&  b) 
{ 
  return ((m_f ? m_f->grad(b) : Box3d(0.0)) + 
          (m_g ? m_g->grad(b) : Box3d(0.0))); 
}

IMatrix3d Sum::hess(const Box<double>&  b) 
{ 
  return ((m_f ? m_f->hess(b) : IMatrix3d(0.0)) + 
          (m_g ? m_g->hess(b) : IMatrix3d(0.0))); 
}

