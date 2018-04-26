/*****
 * @file Plane.cpp
 *
 * This is the implementation of the Plane class. The implicit surface
 * function of a plane is defined as:
 *
 * f(x,y,z) = Ax + By + Cz + D, or
 * f(v) = (v - p) . N
 *
 * where v is the input vector, p is a point on the plane, and N is the
 * unit surface normal. The second method is used here. The gradient is:
 *
 * g(x,y,z) = (Nx, Ny, Nz)
 *
 * where Nx, Ny, and Nz are the unit normal components. The Hessian is a 
 * zero matrix.
 *
 * @author: Ed Bachta
 */

#include "Plane.h"
#include "Surface/gmTNTconvert.h"

REGISTER_IMPLICIT(Plane,"Geometric:Plane");

/**
 * Initialization routine called by constructors.  This function allows all
 * of the main object init code to be in one place.
 */
void Plane::init(const gmVector3& n, double d)
{
  m_n = n;   // Normal vector of plane
  m_n.normalize();         // Make unit length
  m_d = d;   // Shortest distance from origin to plane
}

/**
 * Plane()
 *
 * Creates the default x-y plane
 */
Plane::Plane() :
Geometric()
{
  init(gmVector3(0,0,0),0);
} 

/**
 * Plane (double, double, double, double)
 *
 * Creates the plane described by the general equation (ax + by + cz + d = 0)
 */
Plane::Plane(double a, double b, double c, double d) :
Geometric()
{
  gmVector3 vec = gmVector3(a,b,c);
  init(vec,(-d / vec.length()));
} 

/**
 * Plane (const gmVector3&, double)
 *
 * Creates the plane described by n (normal) and d (shortest distance
 * from the origin to the plane).
 */
Plane::Plane(const gmVector3& n, double d) :
Geometric()
{
  init(n,d);
} 

/**
 * Plane (const gmVector3&, gmVector3)
 *
 * Creates the plane described by n (normal) and p (vector to
 * a point on the plane from the origin)
 */
Plane::Plane(const gmVector3& n, const gmVector3& p) :
Geometric()
{
  m_n = n;   // Normal vector of plane
  m_n.normalize();         // Make unit length
  m_d = dot(p,m_n);   // Shortest distance from origin to plane
 }

/**
 * Convenience function for the calculation of a LRP over a time interval
 * for a particular m_n.  We use the "current" and the "last" coefficients
 * to calculate an interval of the coefficients over a given time interval.
 * This is used for 4D critical point finding.
 */
Box3d Plane::m_nlrp(Intervald t)
{
  double* m_nold = (double *) malloc(qlen() * sizeof(double));
  getqold(m_nold);
  Box3d retval;
  
  // m_n resides in q[0],q[1],q[2] = m_nold[0..2]
  for (int i = 0; i < 3; i++)
    retval[i] = t * (m_n[i] - m_nold[i]) + Intervald(m_nold[i]);

  delete m_nold;
  return retval;
}

Intervald Plane::m_dlrp(Intervald t)
{
  double* m_dold = (double *) malloc(qlen() * sizeof(double));
  getqold(m_dold);
  // m_d resides in q[3] = m_dold[3]
  Intervald retval =  t * (m_d - m_dold[3]) + Intervald(m_d);
  delete m_dold;
  return retval;
}

#ifndef INTERVAL_EVAL_ONLY
double Plane::proc(const gmVector3 & x) 
{
  gmVector3 p = m_n * m_d;

  return dot(x-p, m_n);
} 

gmVector3 Plane::grad(const gmVector3 & x) 
{
  gmVector3 g(m_n[0], m_n[1], m_n[2]);
  
  return g;
}

gmMatrix3 Plane::hess(const gmVector3 & x)
{
  return gmMatrix3();
}
#endif

Intervald Plane::proc(const Box<double>& x)
{
  Box3d x1;
  Box3d p;
  Box3d diff;
  Box3d mn;

  if (x.size() == 3)
    {
      x1 = x;
      mn = Box<double>(convert(m_n));
      p = Box<double>(convert(m_n * m_d));
    }
  else
    {
      x1 = Box3d(x[0],x[1],x[2]);
      mn = m_nlrp(x[3]);
      p = m_nlrp(x[3]) * m_dlrp(x[3]);
    }

  diff = x1 - p;
  return dot(diff,mn);
}

Box3d Plane::grad(const Box<double>& x)
{
  Box3d g;

  if (x.size() == 3)
    g = Box<double>(convert(m_n));
  else
    g = m_nlrp(x[3]);

  return g;
}

IMatrix3d Plane::hess(const Box<double>& x)
{
  return IMatrix3d();
}

/**
 * @todo Fix me - but may not be calculatable.
 */
Intervald Plane::proct(const Box<double>& x)
{
  return Intervald(0.0);
}

/**
 * @todo Fix me - but may not be calculatable.
 */
Box3d Plane::gradt(const Box<double>& x)
{
  return Box3d(0.0);
}

/**
 * procq(const gmVector3&, double*)
 * Returns df/dq
 */
void Plane::procq(const gmVector3 & x, double* dfdq) 
{
  dfdq[0] = x[0] - 2 * m_d * m_n[0];
  dfdq[1] = x[1] - 2 * m_d * m_n[1];
  dfdq[2] = x[2] - 2 * m_d * m_n[2];
  dfdq[3] = - (m_n[0]*m_n[0]) 
            - (m_n[1]*m_n[1]) 
            - (m_n[2]*m_n[2]);
}

/**
 * getq(double*)
 * Loads q with parameters
 */
void Plane::getq(double* q) 
{
  q[0] = m_n[0];
  q[1] = m_n[1];
  q[2] = m_n[2];
  q[3] = m_d;
} 

/**
 * _setq(double*)
 * Loads parameters from q
 */
void Plane::_setq(double* q) 
{
  m_n[0] = q[0]; 
  m_n[1] = q[1]; 
  m_n[2] = q[2];
  m_d = q[3];
}

/**
 * getqname(char**)
 * Returns a list of parameter names
 */
void Plane::getqname(char** qn) 
{
  qn[0] = "Normal X";
  qn[1] = "Normal Y";
  qn[2] = "Normal Z";
  qn[3] = "Displacement";
}


const char ** Plane::getPixmapXPM(const int& size) const
{
  if (size <= 16)
    return (const char **)plane_pixmap16;
  else if (size <= 32)
    return (const char **)plane_pixmap32;
  else
    return (const char **)plane_pixmap48;
}

