/**
 * Implementation of an offset.
 * @file Offset.cpp
 * @date Fall 2000
 * @author Ed Bachta
 */

#include "Offset.h"

REGISTER_IMPLICIT(Offset,"UnaryOp:Offset");

/**
 *  Called by the various constructors to allow for a single location for
 *  the init of a Offset object.
 */
void Offset::init(Implicit* f, double r)
{
  m_f = f;
  m_r = r;
}

/**
 * Default constructor.
 */
Offset::Offset()
{
  init(NULL,0.0);
} 

/**
 * Explicit constructor.
 * @param f Surface to offset.
 * @param r Value of offset.
 */
Offset::Offset(Implicit *f, double r)
{
  init(f,r);
} 

/**
 * Convenience function for the calculation of a LRP over a time interval
 * for a particular m_r.  We use the "current" and the "last" coefficients
 * to calculate an interval of the coefficients over a given time interval.
 * This is used for 4D critical point finding.
 */
Intervald Offset::m_rlrp(Intervald t)
{
  double* m_rold = (double *) malloc(qlen() * sizeof (double));
  getqold(m_rold);
  // m_r resides in q[0] = m_rold[0]
  Intervald retval =  t * (m_r - m_rold[0]) + Intervald(m_rold[0]);
  delete m_rold;
  return retval;
}

/**
 * Convenience function to calculate the difference between the current Q
 * and the old Q.  This is used by proct() and gradt().
 */
double Offset::m_rdiff(void)
{
  double* m_rold = (double *) malloc(qlen() * sizeof (double));
  getqold(m_rold);
  // m_r resides in q[0] = m_rold[0]
  double retval =  m_r - m_rold[0];
  delete m_rold;
  return retval;
}

#ifndef INTERVAL_EVAL_ONLY
/**
 * Evaluation of function.
 * @param   x Point to evaluate at.
 * @returns Function value.
 */
double Offset::proc(const gmVector3 & x)
{
  return (m_f ? m_f->proc(x) - m_r : 0.0);
}

/**
 * Evaluation of gradient.
 * @param   x Point to evaluate at.
 * @returns The gradient at x.
 */
gmVector3 Offset::grad(const gmVector3 & x)
{
  return (m_f ? m_f->grad(x) : gmVector3());
}

/**
 * Evaluation of Hessian.
 * @param   x Point to evaluate at.
 * @returns The Hessian at x.
 */
gmMatrix3 Offset::hess(const gmVector3 & x)
{
  return (m_f ? m_f->hess(x) : gmMatrix3());
}
#endif

Intervald Offset::proc(const Box<double>&  b)
{
  Intervald retval(0.0);
  Intervald m_rtemp;

  if (m_f)
    {
      if (b.size() == 3)
        m_rtemp = m_r;
      else
        m_rtemp = m_rlrp(b[3]);
      retval = m_f->proc(b) - m_rtemp;
    }

  return retval;
}

Box3d Offset::grad(const Box<double>&  b)
{
  return (m_f ? m_f->grad(b) : Box3d(0.0));
}

IMatrix3d Offset::hess(const Box<double>&  b)
{
  return (m_f ? m_f->hess(b) : IMatrix3d(0.0));
}

Intervald Offset::proct(const Box<double>&  b)
{
  Intervald retval(0.0);
  Intervald m_rtemp;

  if (m_f)
    {
      m_rtemp = Intervald(m_rdiff());
      retval = m_f->proct(b) - m_rtemp;
    }

  return retval;
}

Box3d Offset::gradt(const Box<double>&  b)
{
  return (m_f ? m_f->gradt(b) : Box3d(0.0));
}

/**
 * Retreives parameters.
 * @param q Parameter array.
 */
void Offset::getq(double* q)
{
  q[0] = m_r;
  
  // Check for child
  if (m_f)
    m_f->getq(&q[1]);
}

/**
 * Assigns parameters.
 * @param q Parameter array.
 */
void Offset::_setq(double* q)
{
  m_r = q[0];
  
  // Check for child
  if (m_f)
    m_f->_setq(&q[1]);
}

/**
 * Evaluation of dFdq.
 * @param x    Point to evaluate at.
 * @param dfdq Array representing dfdq.
 */
void Offset::procq(const gmVector3 & x, double* dfdq)
{
  if (m_f)
    {
      dfdq[0] = -1;
      m_f->procq(x, &dfdq[1]);
    }
} 

/**
 * Returns the # of parameters.
 */
unsigned int Offset::qlen()
{
  int retval = 1;
  if (m_f)
    retval += m_f->qlen();
  return retval;
}

/**
 * Retreives parameter names.
 * @param qn Paramater name array.
 */
void Offset::getqname(char** qn)
{
  qn[0] = "Offset";

  UnaryOp::getqname(qn);
} 

const char ** Offset::getPixmapXPM(const int& size) const
{
  if (size <= 16)
    return (const char **)offset_pixmap16;
  else if (size <= 32)
    return (const char **)offset_pixmap32;
  else
    return (const char **)offset_pixmap48;
}

