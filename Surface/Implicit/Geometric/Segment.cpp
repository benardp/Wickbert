/*****
 * Implementation of the geometric segment.
 * @file Segment.cpp
 * @date Fall 200
 * @author Ed Bachta
 */

#include "Segment.h"
#include "Surface/gmTNTconvert.h"

// Registration
REGISTER_IMPLICIT(Segment,"Geometric:Segment");

/**
 *  Called by the various constructors to allow for a single location for
 *  the init of a Segment object.
 *  @param a First endpoint.
 *  @param b Second endpoint.
 */
void Segment::init(const gmVector3& a, gmVector3 b)
{
  m_a = a;
  m_b = b;
}

/// Default constructor.
Segment::Segment() 
{
  init(gmVector3(0.0, 0.0, 0.5),gmVector3(0.0, 0.0, -0.5));
} 

/**
 * Explicit constructor.
 * @param a First endpoint.
 * @param b Second endpoint.
 */
Segment::Segment(const gmVector3 &a,const gmVector3& b) 
{
  init (a,b);
}

/**
 * Convenience function for the calculation of a LRP over a time interval
 * for a particular m_a.  We use the "current" and the "last" coefficients
 * to calculate an interval of the coefficients over a given time interval.
 * This is used for 4D critical point finding.
 */
Box3d Segment::m_alrp(Intervald t)
{
  double* m_aold = (double *) malloc(qlen() * sizeof(double));
  getqold(m_aold);
  Box3d retval;

  // m_a resides in q[0],q[1],q[2] = m_aold[0..2]
  for (int i = 0; i < 3; i++)
    retval[i] = t * (m_a[i] - m_aold[i]) + Intervald(m_aold[i]);

  delete m_aold;
  return retval;
}

Box3d Segment::m_blrp(Intervald t)
{
  double* m_bold = (double *) malloc(qlen() * sizeof(double));
  getqold(m_bold);
  Box3d retval;

  // m_b resides in q[3],q[4],q[5] = m_bold[3..5]
  for (int i = 0; i < 3; i++)
    retval[i] = t * (m_b[i] - m_bold[i+3]) + Intervald(m_b[i]);

  delete m_bold;
  return retval;
}

#ifndef INTERVAL_EVAL_ONLY
/**
 * Evaluation of function.
 * @param x Position at which to evaluate.
 * @returns Distance to segment.
 */
double Segment::proc(const gmVector3 & x) 
{
  double retval = 0.0;                 // The value to be returned.
  gmVector3 ab = m_b - m_a;           // The vector from a to b.
  gmVector3 ax = x - m_a;             // The vector from a to x.
  double ablen = ab.length();
  ab.normalize(); // ab is now a simple unit direction vector
  double p = dot(ax,ab);

  if (p <= 0) 
    { // x's projection onto ab extends to or past the a end of the segment 
      retval = ax.length();
    } 
  else if (p >= ablen) 
    { // x's projection onto ab extends to or past the b end of the segment
      retval = (x - m_b).length();
    } 
  else 
    { // x's projection onto ab lies somewhere between a and b
      retval = sqrt(ax.lengthSquared() - p*p);
    }
  return retval;
}
#endif

Intervald Segment::proc(const Box<double>& x)
{
  Intervald retval;
  Box3d x1;
  Box3d ab;
  Box3d ax;
  Box3d m_atemp;
  Box3d m_btemp;
  if (x.size() == 3)
    {
      x1 = x;
      m_atemp = Box<double>(convert(m_a));
      m_btemp = Box<double>(convert(m_b));
      ab = m_btemp - m_atemp;
      ax = x - m_atemp;
    }
  else
    {
      x1 = Box3d(x[0],x[1],x[2]);
      m_atemp = m_alrp(x[3]);
      m_btemp = m_blrp(x[3]);
      ab = m_btemp - m_atemp;
      ax = x - m_atemp;
    }
  Intervald ablen = ab.length();
  ab.normalize();
  Intervald p = dot(ax,ab);

  //  FIX ME RIGHT!!!  This is just a temporary hack right now. 
  if (p.center() <= 0)
    { // x's projection onto ab extends to or past the a end of the segment 
      retval = ax.length();
    }
  else if (p.center() >= ablen)
    { // x's projection onto ab extends to or past the b end of the segment
      Box3d tempor = x1 - m_btemp;
      retval = tempor.length();
    }
  else
    { // x's projection onto ab lies somewhere between a and b
      retval = (ax.lengthSquared() - p*p).sqrt();
    }

  return retval;
}

#ifndef INTERVAL_EVAL_ONLY
/**
 * Evaluation of gradient
 * @param   x Position at which to evaluate.
 * @returns Gradient vector.
 */
gmVector3 Segment::grad(const gmVector3 & x) 
{
  gmVector3 g;
  gmVector3 ab = m_b - m_a;           // vector from a to b
  gmVector3 ax = x - m_a;             // vector from a to x
  double ablen = ab.length();
  ab.normalize(); // ab is now a simple unit direction vector
  double p = dot(ax, ab);
  
  if (p <= 0) 
    { // x's projection onto ab extends to or past the a end of the segment
      g = ax;
    } 
  else if (p >= ablen) 
    { // x's projection onto ab extends to or past the b end of the segment
      g = x - m_b;
    } 
  else 
    { // x's projection onto ab lies somewhere between a and b
      g = ax - (ab*p); // find normal vector
    }  

  return g.normalize();
}
#endif

Box3d Segment::grad(const Box<double>& x)
{
  Box3d g;
  Box3d x1;
  Box3d ab;
  Box3d ax;
  Box3d m_atemp;
  Box3d m_btemp;
  if (x.size() == 3)
    {
      x1 = x;
      m_atemp = Box<double>(convert(m_a));
      m_btemp = Box<double>(convert(m_b));
      ab = m_btemp - m_atemp;
      ax = x - m_atemp;
    }
  else
    {
      x1 = Box3d(x[0],x[1],x[2]);
      m_atemp = m_alrp(x[3]);
      m_btemp = m_blrp(x[3]);
      ab = m_btemp - m_atemp;
      ax = x - m_atemp;
    }
  Intervald ablen = ab.length();
  ab.normalize();
  Intervald p = dot(ax,ab);

  //  FIX ME RIGHT!!!  This is just a temporary hack right now. 
  if (p.center() <= 0) 
    { // x's projection onto ab extends to or past the a end of the segment
      g = ax;
    } 
  else if (p.center() >= ablen) 
    { // x's projection onto ab extends to or past the b end of the segment
      g = x1 - m_btemp;
    } 
  else 
    { // x's projection onto ab lies somewhere between a and b
      g = ax - (ab*p); // find normal vector
    }  

  return g.normalize();
}

#ifndef INTERVAL_EVAL_ONLY
gmMatrix3 Segment::hess(const gmVector3 & x)
{
  // Need to fix this!
  gmMatrix3 h; // Initialized to zero

  return h;
}
#endif

IMatrix3d Segment::hess(const Box<double>& x)
{
  // Need to fix this!
  IMatrix3d h = IMatrix3d(0.0);

  return h;
}

/**
 * @todo Fix me - but may not be calculatable.
 */
Intervald Segment::proct(const Box<double>& x)
{
  return Intervald(0.0);
}

/**
 * @todo Fix me - but may not be calculatable.
 */
Box3d Segment::gradt(const Box<double>& x)
{
  return Box3d(0.0);
}

/**
 * getq(double*)
 * Loads parameters into q
 */
void Segment::getq(double* q) 
{
  q[0] = m_a[0];
  q[1] = m_a[1];
  q[2] = m_a[2];
  q[3] = m_b[0];
  q[4] = m_b[1];
  q[5] = m_b[2];
}

/**
 * _setq(double*)
 * Sets attributes from q
 */
void Segment::_setq(double* q) 
{
  m_a[0] = q[0]; 
  m_a[1] = q[1]; 
  m_a[2] = q[2];
  m_b[0] = q[3]; 
  m_b[1] = q[4]; 
  m_b[2] = q[5];
}


/**
 * getqname(char**)
 * Returns a list of parameter names
 */
void Segment::getqname(char** qn) 
{
  qn[0] = "X0 coord";
  qn[1] = "Y0 coord";
  qn[2] = "Z0 coord";
  qn[3] = "X1 coord";
  qn[4] = "Y1 coord";
  qn[5] = "Z1 coord";
}

