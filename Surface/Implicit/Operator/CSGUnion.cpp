//*****************************************************************************
//	CSGUnion.h: Implementation file for the CSGUnion class.
//	Author: Xinlai Ni
//	1/13/2006
//*****************************************************************************

#include <string.h>
#include "CSGUnion.h"

REGISTER_IMPLICIT(CSGUnion, "BinaryOp:CSGUnion");


void CSGUnion::procq(const gmVector3 & x, double* q)
{
	if( m_f != NULL && m_g != NULL ) {
		m_f->procq(x, q);
		m_g->procq(x, q	+ m_f->qlen());

		double fx = m_f->proc(x);
		double gx = m_g->proc(x);

		double dhdf = hf(fx, gx);
		unsigned int i;
		for(i = 0; i < m_f->qlen(); i++)
			q[i] *= dhdf;
		double dhdg = hg(fx, gx);
		for(; i < m_f->qlen() + m_g->qlen(); i++)
			q[i] *= dhdg;
	}
}

#ifndef INTERVAL_EVAL_ONLY
double CSGUnion::proc(const gmVector3 & x)
{
  if ((m_f!=NULL) && (m_g!=NULL))
    return h(m_f->proc(x),m_g->proc(x));
  else
    return 0.0;
}

gmVector3 CSGUnion::grad(const gmVector3 & x)
{
  if ((m_f!=NULL) && (m_g!=NULL))
    {
      double fx = m_f->proc(x); 
      double gx = m_g->proc(x);
      return hf(fx,gx) * m_f->grad(x) + hg(fx,gx) * m_g->grad(x);
    }
  else
    return gmVector3();
}

gmMatrix3 CSGUnion::hess(const gmVector3 & x)
{
  if ((m_f!=NULL) && (m_g!=NULL))
    {
      double fx = m_f->proc(x); 
      double gx = m_g->proc(x);
      gmVector3 dfx = m_f->grad(x);
      gmVector3 dgx = m_g->grad(x);

      return hff(fx,gx) * outer(dfx,dfx) + 
             hfg(fx,gx) * outer(dfx,dgx) +
             hfg(fx,gx) * outer(dgx,dfx) + 
             hgg(fx,gx) * outer(dgx,dgx) +
             hf(fx,gx)  * m_f->hess(x)   + 
             hg(fx,gx)  * m_g->hess(x);
    }
  else
    return gmMatrix3();
}
#endif