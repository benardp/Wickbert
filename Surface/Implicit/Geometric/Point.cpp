/**
 * Implementation of the geometric point
 * @file Point.cpp
 * @date Fall 2000
 * @author Ed Bachta
 */

#include "Point.h"
#include "Surface/gmTNTconvert.h"

REGISTER_IMPLICIT(Point,"Geometric:Point");

/// Default constructor.
Point::Point() : 
Sphere(0.0) 
{
}

/**
 * Explicit contstructor.
 * @param x Position.
 */
Point::Point(const gmVector3 & x) : 
Sphere(x,0.0) 
{
}

/**
 * Gets parameters.
 * @param q An array for parameters.
 */
void Point::getq(double* q) 
{
  q[0] = m_x[0];
  q[1] = m_x[1];
  q[2] = m_x[2];
}

/**
 * Sets parameters.
 * @param q An array of parameters.
 */
void Point::_setq(double* q) 
{
  m_x[0] = q[0]; 
  m_x[1] = q[1]; 
  m_x[2] = q[2];
} 

/**
 * Get parameter names.
 * @param qn An array for names.
 */
void Point::getqname(char** qn) 
{
  qn[0] = "X coord";
  qn[1] = "Y coord";
  qn[2] = "Z coord";
}

