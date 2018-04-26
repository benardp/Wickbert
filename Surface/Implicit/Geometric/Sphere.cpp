/**
 * Implementation of the geometric sphere.
 * @file Geometric/Sphere.cpp
 * @date Fall 2000
 * @author Ed Bachta
 * @note Changed to geoetric sphere (from algebraic), July 2001, jch
 */

#include "Sphere.h"
#include "Surface/gmTNTconvert.h"

REGISTER_IMPLICIT(Sphere,"Geometric:Sphere");

/**
 *  Called by the various constructors to allow for a single location for
 *  the init of a Sphere object.
 */
void Sphere::init(const gmVector3 & x, double r)
{
  m_x = x;
  m_r = r;
}

/**
 * Creates the default sphere of radius 1 at the origin.
 */
Sphere::Sphere() :
Geometric()
{
  init(gmVector3(0.0,0.0,0.0),1.0);
}

/**
 * Creates a sphere with radius r at the origin.
 */
Sphere::Sphere(double r) :
Geometric()
{
  init(gmVector3(0.0,0.0,0.0),r);
} 

/**
 * Creates a sphere with radius r at position x.
 */
Sphere::Sphere(const gmVector3 & x, double r):
Geometric()
{
  init(x,r);
}

/**
 * Convenience function for the calculation of a LRP over a time interval
 * for a particular m_x.  We use the "current" and the "last" coefficients
 * to calculate an interval of the coefficients over a given time interval.
 * This is used for 4D critical point finding.
 */
Box3d Sphere::m_xlrp(Intervald t)
{
  double* m_xold = (double *) malloc(qlen() * sizeof(double));
  Box3d retval;
  
  // m_x resides in q[0],q[1],q[2] = m_xold[0..2]
  for (int i = 0; i < 3; i++)
    retval[i] = t * (m_x[i] - m_xold[i]) + Intervald(m_xold[i]);

  delete m_xold;
  return retval;
}

Intervald Sphere::m_rlrp(Intervald t)
{
  double* m_rold = (double *) malloc(qlen() * sizeof(double));
  getqold(m_rold);
  // m_r resides in q[3] = m_dold[3]
  Intervald retval = t * (m_r - m_rold[3]) + Intervald(m_r);
  delete m_rold;
  return retval;
}

#ifndef INTERVAL_EVAL_ONLY
/**
 * Evaluates the implicit function of a sphere: 
 * f(x) = ||x - m_x|| - R
 * @param   x Point at which to evaluate function.
 * @returns Function value.
 */
double Sphere::proc(const gmVector3 & x)
{
  return (x - m_x).length() - m_r;
}

/**
 * Evaluates the gradient of a sphere.
 * @param   x Point at which to evaluate gradient.
 * @returns Gradient vector.
 */
gmVector3 Sphere::grad(const gmVector3 & x) 
{
  gmVector3 g = x - m_x;
  double l = g.lengthSquared();
  if (gmIsZero(l)) 
    return gmVector3(0.0,0.0,0.0);
  else 
    return (g / sqrt(l));
} 

/**
 * Returns the Hessian of a Geometric sphere:
 *
 * f(X) = (x*x + y*y + z*z)^1/2 - r
 * df/dx = 2x * 1/2(x*x + y*y + z*z)^-1/2
 * grad f = x/(X.X)^1/2
 * d^2f/dx^2 = 1/(X.X)^1/2 - x^2/(X.X)^3/2
 * d^2f/dxdy = -xy/(X.X)^3/2
 *
 * @param   x Position irrelevant for a sphere.
 * @returns Hessian matrix.
 */
gmMatrix3 Sphere::hess(const gmVector3 & x) 
{
  gmMatrix3 h;  // initialized to zero
  double l = x.lengthSquared();

  if (!gmIsZero(l))
    {
      double li = 1.0/sqrt(l);
      double li3 = li*li*li;

      h[0][0] = l - x[0]*x[0]*li3;
      h[1][1] = l - x[1]*x[1]*li3;
      h[2][2] = l - x[2]*x[2]*li3;
      h[0][1] = h[1][0] = -x[0]*x[1]*li3;
      h[0][2] = h[2][0] = -x[0]*x[2]*li3;
      h[1][2] = h[2][1] = -x[1]*x[2]*li3;
    }

  return h;
} 
#endif

Intervald Sphere::proc(const Box<double>& x) 
{
  Box3d x1;
  Box3d center;
  Box3d diff;
  Intervald m_rtemp;

  if (x.size() == 3)
    {
      x1 = x;
      center = Box<double>(convert(m_x));
      m_rtemp = Intervald(m_r);
    }
  else
    {
      x1 = Box3d(x[0],x[1],x[2]);
      center =  m_xlrp(x[3]);
      m_rtemp = m_rlrp(x[3]);
    }

  diff = x1 - center;
  return diff.length() - m_rtemp;
} 

Box3d Sphere::grad(const Box<double>& x)
{
  Box3d x1;
  Box3d center;
  Box3d g;

  if (x.size() == 3)
    {
      x1 = x;
      center = Box<double>(convert(m_x));
    }
  else
    {
      x1 = Box3d(x[0],x[1],x[2]);
      center = m_xlrp(x[3]);
    }

  g = x1 - center;
  Intervald l = g.lengthSquared();
  if (l.isZero())
    return Box3d(0.0);
  else
    return g / l.sqrt();
}

IMatrix3d Sphere::hess(const Box<double>& x)
{
  IMatrix3d h = IMatrix3d(0.0);
  Box3d x1 = Box3d(x[0],x[1],x[2]);
  Intervald l = x1.lengthSquared();

  if (!l.isZero())  
    {
      Intervald li  = Intervald(1.0) / l.sqrt();
      Intervald li3 = li.pow(3);

      h[0][0] = l - x[0].squared() * li3;
      h[1][1] = l - x[1].squared() * li3;
      h[2][2] = l - x[2].squared() * li3;
      h[0][1] = h[1][0] = -x[0] * x[1] * li3;
      h[0][2] = h[2][0] = -x[0] * x[2] * li3;
      h[1][2] = h[2][1] = -x[1] * x[2] * li3;
    }

  return h;
}

/**
 * @todo Fix me!
 */
Intervald Sphere::proct(const Box<double>& x)
{
  return Intervald(0.0);
}

/**
 * @todo Fix me!
 */
Box3d Sphere::gradt(const Box<double>& x)
{
  return Box3d(0.0);
}

/**
 * Evaluates df/dq
 * f = (x-c).(x-c)^1/2 - r
 * df/dc = -(x-c)/(((x-c).(x-c))^1/2)
 * df/dr = -1
 */
void Sphere::procq(const gmVector3 & x, double* dfdq) 
{
  gmVector3 x_temp = x;
  x_temp -= m_x;

  double l = x_temp.lengthSquared();

  gmVector3 dfdc = -x_temp/sqrt(l);  // what if sqrt(l) near zero?

  dfdq[0] = dfdc[0];
  dfdq[1] = dfdc[1];
  dfdq[2] = dfdc[2];
  dfdq[3] = -1.0;
} 

/**
 * Loads parameters into q.
 */
void Sphere::getq(double* q) 
{
  q[0] = m_x[0];
  q[1] = m_x[1];
  q[2] = m_x[2];
  q[3] = m_r;
}

/**
 * Loads parameters from q.
 */
void Sphere::_setq(double* q) 
{
  m_x[0] = q[0]; 
  m_x[1] = q[1]; 
  m_x[2] = q[2];
  m_r = q[3];
}

/**
 * Returns a list of parameter names.
 * \param qn an array of size qlen of strings listing names of parameters
 */
void Sphere::getqname(char** qn) 
{
  qn[0] = "Center X";
  qn[1] = "Center Y";
  qn[2] = "Center Z";
  qn[3] = "Radius";
}

const char ** Sphere::getPixmapXPM(const int& size) const
{
  if (size <= 16)
    return (const char **)sphere_pixmap16;
  else if (size <= 32)
    return (const char **)sphere_pixmap32;
  else
    return (const char **)sphere_pixmap48;
}

