/**
 * @File Wyvill.cpp
 * Implementation of Wyvill class
 */

#include <math.h>
#include "Wyvill.h"

REGISTER_IMPLICIT(Wyvill,"UnaryOp:Blob:Wyvill");

/**
 *  Called by the various constructors to allow for a single location for
 *  the init of a Wyvill object.
 */
void Wyvill::init(Quadric *f, double r)
{
  m_f = f;
  m_r = r;
}

/// Default constructor.
Wyvill::Wyvill()
{
  init(NULL,1.0);
}

/**
 * Class Constructor:
 * Takes a pointer to an arbitrary quadric 
 * as an argument and initializes Radius to 1
 */
Wyvill::Wyvill(Quadric *f)
{
  init(f,1.0);
}

/**
 * Class Constructor:
 * Takes a pointer to an arbitrary quadric 
 * and a scalar radius r as agruments
 */
Wyvill::Wyvill(Quadric *f, double r)
{
  init(f,r);
}

/**
 * Returns the evaluation of the Polynomial Blob  
 * at the specified 3-space point
 *
 * The Poly Blob f(v) is defined as:
 *    f(v)      = -4/9(r^2/R^2)^3 + 17/9(r^2/R^2)^2 - 22/9(r^2/R^2) + 1
 *    where r^2 = quadric radius function
 */
double Wyvill::kernel(double r2)
{
  double blob;
  double R2 = m_r * m_r;

  // definiton of Wyvill Blob:
  // 
  blob = (-0.44444f / (R2 * R2 * R2)) * (r2 * r2 * r2) + 
         ( 1.88888f / (R2 * R2))      * (r2 * r2)      + 
         (-2.44444f /  R2)            *  r2            + 1.0f;

  return blob;
}

/*****
 * Kernel derivative of the Polynomial Blob is
 * defined as follows:
 *
 *    d/dr^2    = (-12/(9R^6)) ((r^2)^2) + (34/(9R^4)) (r^2) + (-22/(9R^2))  
 *    where r^2 = quadric radius function
 *            R = rest radius
 */
double Wyvill::dkernel(double r2)
{
  double R2 = m_r * m_r;
  double dblob;

  dblob = (-1.33333f / (R2 * R2 *R2)) * (r2 * r2) +
          ( 3.77778f / (R2 * R2)    ) *  r2       +
          (-2.44444f /  R2); 

  return dblob;
}

/**
 * Kernel second derivative of the Polynomial Blob  
 * is defined as follows:
 *
 *    d/d2r^2   = (-24/(9R^6))r^2 + (34/(9R^4))
 *    where r^2 = quadric radius function
 *            R = rest radius
 */
double Wyvill::d2kernel(double r2)
{
  double R2 = m_r * m_r;
  double dblob;

  dblob = (-2.66667f / (R2 * R2 *R2)) *  r2 +
          ( 3.77778f / (R2 * R2)    ); 

  return dblob;
}

Intervald Wyvill::kernel(Intervald r2)
{
  Intervald blob;
  double R2 = m_r * m_r;

  blob = Intervald(-0.44444f / (R2 * R2 * R2)) * (r2 * r2 * r2) + 
         Intervald( 1.88888f / (R2 * R2))      * (r2 * r2)      + 
         Intervald(-2.44444f /  R2)            *  r2            + 
         Intervald(1.0);

  return blob;
}

Intervald Wyvill::dkernel(Intervald r2)
{
  Intervald dblob;
  double R2 = m_r * m_r;

  dblob = Intervald(-1.33333f / (R2 * R2 *R2)) * (r2 * r2) +
          Intervald( 3.77778f / (R2 * R2)    ) *  r2       +
          Intervald(-2.44444f /  R2); 

  return dblob;
}

Intervald Wyvill::d2kernel(Intervald r2)
{
  Intervald dblob;
  double R2 = m_r * m_r;

  dblob = Intervald(-2.66667f / (R2 * R2 *R2)) *  r2 +
          Intervald( 3.77778f / (R2 * R2)    ); 

  return dblob;
}

/**
 * Function qlen returns the parameter list 
 * for the definition of a Blob
 */
unsigned int Wyvill::qlen(void)
{
  // queue length equal to length of 
  // queue for the radius function plus
  // rest radius variable
  int retval = 1;
  if (m_f)
    retval += m_f->qlen();
  return retval;
}

/**
 * Function setq sets the parameters
 * defining the Blob to those specified
 * by array q
 * @note Does not alter q
 */
void Wyvill::_setq(double *q)
{
  m_r = q[0];

  if (m_f == NULL)
    return;

  double *m_fq = new double[m_f->qlen ()];

  // alter queue parameters of the
  // radius function
  for (unsigned int i = 1; i < qlen (); i++)
    m_fq[i - 1] = q[i];
  m_f->_setq(m_fq);

  delete [] m_fq;
}

/**
 * Function getq returns an array q of the 
 * parameter values defining the present
 * Blob
 *
 * Array length required to be at least 
 * this->qlen () in size
 */
void Wyvill::getq(double *q)
{
  q[0] = m_r;

  if (m_f == NULL)
    return;

  double *m_fq = new double[m_f->qlen ()];
  m_f->getq (m_fq);

  // add queue parameters specified
  // by the radius function
  for (unsigned int i = 1; i < qlen (); i++)
    q[i] = m_fq[i - 1];

  delete [] m_fq;
}

/**
 * Function getqname returns an array of 
 * strings defining in english the names
 * of parameters for the Blob
 *
 * Array required to have at least qlen ()
 * pointers to strings
 */
void Wyvill::getqname(char** qn)
{
  qn[0] = "Rest radius R";

  if (m_f != NULL)
    {
      char **m_fqname = new char*[m_f->qlen ()];
      m_f->getqname(m_fqname);

      for (unsigned int i = 1; i < qlen (); i++)
        qn[i] = m_fqname[i - 1];

	  //MEMORY LEAK ???? I feel like this is also causing memory problems!
	  //getqname creates the entries (RBF) other getqname functions do not! (compactRBF, Plane)
	  //I feel like there is a major problem here...
	  //Anyway we would have to delete all pointers, then the pointer m_fqname -- EE
      delete [] m_fqname;
    }
}

/**
 * Derivatives of parameters defined 
 * as follows:
 *
 *    d/dR   = (24/(9R^7)) ((r^2)^3) - (68/(9R^5)) ((r^2)^2) + (44/(9R^3))r^2  
 *    d/dr^2 = (-12/(9R^6)) ((r^2)^2) (df/dq r^2) + 
 *             ( 34/(9R^4))  (r^2)    (df/dq r^2) + 
 *             (-22/(9R^2))           (df/dq r^2)
 *    where r^2 = quadric radius function
 *            R = rest radius
 */
void Wyvill::procq(const gmVector3 & x, double* dfdq)
{
  if (m_f != NULL)
  {
    double R2 = m_r * m_r;
    double r2 = m_f->proc (x);

    dfdq[0] = ( 2.66667f / (R2 * R2 *R2 * m_r)) * (r2 * r2 * r2) +
              (-7.55556f / (R2 * R2     * m_r)) * (r2 * r2     ) +
              ( 4.88889f / (R2          * m_r)) *  r2; 

    double *m_fdfdq = new double[m_f->qlen()];
    m_f->procq (x, m_fdfdq);

    for (unsigned int i = 1; i < qlen (); i++)
      dfdq[i] = (-1.33333f / (R2 * R2 *R2)) * (r2 * r2) * m_fdfdq[i - 1] +
                ( 3.77778f / (R2 * R2)    ) *  r2       * m_fdfdq[i - 1] +
                (-2.44444f /  R2)                       * m_fdfdq[i - 1]; 

    delete [] m_fdfdq;
  }
}

