#include "Blend.h"

/**
 *  Called by the various constructors to allow for a single location for
 *  the init of a Sphere object.
 */
void Blend::init(Implicit* f, Implicit* g)
{
  m_f = f;
  m_g = g;
  m_r1 = 1.0;
  m_r2 = 1.0;
}

#ifndef INTERVAL_EVAL_ONLY
double Blend::proc(const gmVector3 & x)
{
  if ((m_f!=NULL) && (m_g!=NULL))
    return h(m_f->proc(x),m_g->proc(x));
  else
    return 0.0;
}

gmVector3 Blend::grad(const gmVector3 & x)
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

gmMatrix3 Blend::hess(const gmVector3 & x)
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

Intervald Blend::proc(const Box<double>& x)
{
  if ((m_f!=NULL) && (m_g!=NULL))
    return h(m_f->proc(x),m_g->proc(x));
  else
    return Intervald(0.0);
}

Box3d Blend::grad(const Box<double>& x)
{
  if ((m_f!=NULL) && (m_g!=NULL))
    {
      Intervald fx = m_f->proc(x);
      Intervald gx = m_g->proc(x);
      return hf(fx,gx) * m_f->grad(x) + hg(fx,gx) * m_g->grad(x);
    }
  else
    return Box3d(0.0);
}

IMatrix3d Blend::hess(const Box<double>& x)
{
  if ((m_f!=NULL) && (m_g!=NULL))
    {
      Intervald fx = m_f->proc(x); 
      Intervald gx = m_g->proc(x);
      Box3d dfx = m_f->grad(x);
      Box3d dgx = m_g->grad(x);

      return hff(fx,gx) * outer(dfx,dfx) + 
             hfg(fx,gx) * outer(dfx,dgx) +
             hfg(fx,gx) * outer(dgx,dfx) + 
             hgg(fx,gx) * outer(dgx,dgx) +
             hf(fx,gx)  * m_f->hess(x)   + 
             hg(fx,gx)  * m_g->hess(x);
    }
  else
    return IMatrix3d(0.0);
}

void Blend::_setq(double* q)
{
  m_r1 = q[0];
  m_r2 = q[1];

  if ((m_f!=NULL) && (m_g!=NULL))
    {
      m_f->_setq(&q[2]);
      m_g->_setq(&q[2+m_f->qlen()]);
    } 
}

void Blend::getq(double* q)
{  
  q[0] = m_r1;
  q[1] = m_r2;

  if ((m_f!=NULL) && (m_g!=NULL))
    {
      m_f->getq(&q[2]);
      m_g->getq(&q[2+m_f->qlen()]);
    }  
}

void Blend::procq(const gmVector3 & x, double* q)
{
  if ((m_f!=NULL) && (m_g!=NULL))
    {
      q[0] = hr1(m_f->proc(x),m_g->proc(x));
      q[1] = hr2(m_f->proc(x),m_g->proc(x));

      m_f->procq(x, &q[2]);
      m_g->procq(x, &q[2+m_f->qlen()]);
    }
}

void Blend::getqname(char** name)
{
  name[0] = "Radius 1";
  name[1] = "Radius 2";

  if ((m_f!=NULL) && (m_g!=NULL))
    {
      m_f->getqname(&name[2]);
      m_g->getqname(&name[2+m_f->qlen()]);
    }
}

unsigned int Blend::qlen(void)
{
  int retval = 2;

  if ((m_f!=NULL) && (m_g!=NULL))
    retval += m_f->qlen() + m_g->qlen();

  return retval;
}

