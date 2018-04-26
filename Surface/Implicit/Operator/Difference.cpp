/**
 * @file Difference.cpp
 * Implemenation of difference operation
 */
#include "Difference.h"

REGISTER_IMPLICIT(Difference,"BinaryOp:Difference");

#ifndef INTERVAL_EVAL_ONLY
double Difference::proc(const gmVector3 & x) 
{
  return ((m_f ? m_f->proc(x) : 0.0) - 
          (m_g ? m_g->proc(x) : 0.0));
}

gmVector3 Difference::grad(const gmVector3 & x) 
{
  return ((m_f ? m_f->grad(x) : gmVector3()) - 
          (m_g ? m_g->grad(x) : gmVector3()));
}

gmMatrix3 Difference::hess(const gmVector3 & x) 
{
  return ((m_f ? m_f->hess(x) : gmMatrix3()) - 
          (m_g ? m_g->hess(x) : gmMatrix3()));
}
#endif

Intervald Difference::proc(const Box<double>&  b) 
{ 
  return ((m_f ? m_f->proc(b) : Intervald(0.0)) - 
          (m_g ? m_g->proc(b) : Intervald(0.0)));
}

Box3d Difference::grad(const Box<double>&  b) 
{
  return ((m_f ? m_f->grad(b) : Box3d(0.0)) - 
          (m_g ? m_g->grad(b) : Box3d(0.0))); 
}

IMatrix3d Difference::hess(const Box<double>&  b) 
{
  return ((m_f ? m_f->hess(b) : IMatrix3d(0.0)) - 
          (m_g ? m_g->hess(b) : IMatrix3d(0.0))); 
}

