/**
 * @file Blinn.cpp
 * Implementation of Blinn class
 */

#include <math.h>
#include "Blinn.h"

REGISTER_IMPLICIT(Blinn,"UnaryOp:Blob:Blinn");

/**
 *  Called by the various constructors to allow for a single location for
 *  the init of a Blinn object.
 */
void Blinn::init(Quadric *f, double b, double r)
{
  m_f = f;
  m_r = r;
  m_b = b;
}

/// Default constructor.
Blinn::Blinn()
{
  init(NULL,1.0,1.0);
}

/**
 * Class Constructor
 * Takes a pointer to an arbitrary quadric 
 * as an argument and initializes Blobbiness 
 * and Radius to 1
 */
Blinn::Blinn(Quadric *f)
{
  init(f,1.0,1.0);
}

/**
 * Class Constructor:
 * Takes a pointer to an arbitrary quadric 
 * and scalar Blobbiness b as agruments and
 * initializes Radius to 1
 */
Blinn::Blinn(Quadric *f, double b)
{
  init(f,b,1.0);
}

/**
 * Class Constructor:
 * Takes a pointer to an arbitrary quadric, 
 * a scalar Blobbiness b, and a scalar Blob
 * Radius r as agruments
 */
Blinn::Blinn(Quadric *f, double b, double r)
{
  init(f,b,r);
}

/**
 * Returns the evaluation of the Blob at 
 * the specified 3-space point
 *
 * The Blinn Blob f(v) is defined as:
 *    f(v) = exp (-(B/R^2)r^2 + B)
 *
 *    where r^2 = quadric f(v)
 */
double Blinn::kernel(double r2)
{
  double blob;

  // definiton of Blinn Blob:
  blob = exp ((-1.0f * (m_b / (m_r * m_r)) * r2 + m_b));        

  return blob;
}

/**
 * Kernel derivative of the Blinn Blob is
 * defined as follows:
 *
 *    d/dr^2 = (-B/R^2) f(r^2)
 *
 *    where f(r^2) = kernel eval'd at r^2
 *           B = blobbiness
 *           R = rest radius
 */
double Blinn::dkernel(double r2)
{
  double c = -1.0f * m_b / (m_r * m_r); 
  return c * kernel(r2);
}

/**
 * Kernel second derivative of the Blinn Blob  
 * is defined as follows:
 *
 *    d/d2r^2 = (-B/R^2)(-B/R^2) f(r^2)
 *
 *    where f(r^2) = kernel eval'd at r^2
 */
double Blinn::d2kernel(double r2)
{
  double c = -1.0f * m_b / (m_r * m_r); 
  return c * dkernel(r2);
}

Intervald Blinn::kernel(Intervald r2)
{
  Intervald blob;

  blob = Intervald(-1.0) * Intervald(m_b / (m_r * m_r)) * r2 + Intervald(m_b);
  return blob.exp();
}

Intervald Blinn::dkernel(Intervald r2)
{
  double c = -1.0f * m_b / (m_r * m_r);
  return kernel(r2) * Intervald(c);
}

Intervald Blinn::d2kernel(Intervald r2)
{
  double c = -1.0f * m_b / (m_r * m_r);
  return dkernel(r2) * Intervald(c);
}

/**
 * Function qlen returns the parameter list 
 * for the definition of a Blob
 */
unsigned int Blinn::qlen(void)
{
  // queue length equal to length of 
  // queue for the radius function plus
  // rest radius and blobbiness variables
  int retval = 2;
  if (m_f)
    retval += m_f->qlen();
  return retval;
}

/**
 * Function _setq sets the parameters
 * defining the Blob to those specified
 * by array q
 * @note Does not alter q
 */
void Blinn::_setq(double *q)
{
  m_r = q[0];
  m_b = q[1];

  if (m_f == NULL)
    return;

  double *m_fq = new double[m_f->qlen()];

  // alter queue parameters of the
  // radius function
  for (unsigned int i = 2; i < qlen (); i++)
    m_fq[i - 2] = q[i];

  m_f->_setq(m_fq);

  delete[] m_fq;
}

/**
 * Function getq returns an array q of the 
 * parameter values defining the present
 * Blob
 *
 * Array length required to be at least 
 * this->qlen () in size
 */
void Blinn::getq(double *q)
{
  q[0] = m_r;
  q[1] = m_b;

  if (m_f == NULL)
    return;

  double *m_fq = new double[m_f->qlen()];
  m_f->getq(m_fq);

  // add queue parameters specified
  // by the radius function
  for (unsigned int i = 2; i < qlen (); i++)
    q[i] = m_fq[i - 2];

  delete[] m_fq;
}

/**
 * Function getqname returns an array of 
 * strings defining in english the names
 * of parameters for the Blob
 *
 * Array required to have at least qlen ()
 * pointers to strings
 */
void Blinn::getqname(char** qn)
{
  qn[0] = "Rest radius R";
  qn[1] = "Blobbiness  B"; 

  if (m_f != NULL)
    {
      char **m_fqname = new char*[m_f->qlen()];
      m_f->getqname(m_fqname);

      for (unsigned int i = 2; i < qlen (); i++)
        qn[i] = m_fqname[i - 2];

      delete[] m_fqname;
   }
}

/**
 * Derivatives of parameters defined 
 * as follows:
 *
 *    d/dB   = f(v)(-(1/R^2)r^2 + 1) 
 *    d/dR   = f(v)((2B/R^3)r^2)
 *    d/dr^2 = f(v)(-(B/R^2)(df/dq r^2))
 *
 *    where f(v) = exp (-(B/R^2)r^2 + B)
 *           r^2 = quadric radius function
 */
void  Blinn::procq(const gmVector3 & x, double* dfdq)
{
  // an operand surface is defined
  if (m_f != NULL)
    {
      double R2   = m_r * m_r;
      double r2   = m_f->proc(x);
      double blob = m_f->proc(x);

      dfdq[0] = blob * (( 2.0f * r2 * m_b) / (R2 * m_r));
      dfdq[1] = blob * ((-1.0f * r2)                 /  R2+ 1.0f);

      double *m_fdfdq = new double[m_f->qlen ()];
      m_f->procq (x, m_fdfdq);

      for (unsigned int i = 2; i < qlen (); i++)
        dfdq[i] = blob * ((-1.0f * m_fdfdq[i - 2] * m_b) / R2);

      delete[] m_fdfdq;
    }
}

