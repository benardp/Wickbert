#include "Quadric.h"

// Add Quadric to the ImplicitFactory
REGISTER_IMPLICIT(Quadric,"Algebraic:Quadric");

/// Constructs a x^2 + e y^2 + h z^2 + j.
Quadric::Quadric(double a, double e, double h, double j) : Algebraic(2)
{
  setCoef(2,0,0,a);
  setCoef(0,2,0,e);
  setCoef(0,0,2,h);
  setCoef(0,0,0,j);
}

/** Constructs a x^2 + 2b xy + 2c xz + 2d x + e y^2 + 
 *             2f yz + 2g y + h z^2 + 2i z + j
 */
Quadric::Quadric(double a, double b, double c, double d, double e, double f, 
                 double g, double h, double i, double j) : Algebraic(2)
{
  setCoef(2,0,0,a);
  setCoef(1,1,0,2.0*b);
  setCoef(1,0,1,2.0*c);
  setCoef(1,0,0,2.0*d);
  setCoef(0,2,0,e);
  setCoef(0,1,1,2.0*f);
  setCoef(0,1,0,2.0*g);
  setCoef(0,0,2,h);
  setCoef(0,0,1,2.0*i);
  setCoef(0,0,0,j);
}

Quadric::Quadric(double q[10]) : Algebraic(2)
{
  setCoef(2,0,0,q[0]);
  setCoef(1,1,0,2.0*q[1]);
  setCoef(1,0,1,2.0*q[2]);
  setCoef(1,0,0,2.0*q[3]);
  setCoef(0,2,0,q[4]);
  setCoef(0,1,1,2.0*q[5]);
  setCoef(0,1,0,2.0*q[6]);
  setCoef(0,0,2,q[7]);
  setCoef(0,0,1,2.0*q[8]);
  setCoef(0,0,0,q[9]);
}

gmMatrix4 Quadric::Q()
{
  gmMatrix4 qMatrix(getCoef(2,0,0), 0.5*getCoef(1,1,0), 
                    0.5*getCoef(1,0,1), 0.5*getCoef(1,0,0),
                    0.5*getCoef(1,1,0), getCoef(0,2,0), 
                    0.5*getCoef(0,1,1), 0.5*getCoef(0,1,0),
                    0.5*getCoef(1,0,1), 0.5*getCoef(0,1,1), 
                    getCoef(0,0,2), 0.5*getCoef(0,0,1),
                    0.5*getCoef(1,0,0), 0.5*getCoef(0,1,0), 
                    0.5*getCoef(0,0,1), getCoef(0,0,0));
  return qMatrix;
}

void Quadric::Transform(gmMatrix4 t)
{
  gmMatrix4 tAdjoint, qTransform;

  tAdjoint = t.adjoint();
  qTransform = Q() * tAdjoint.transpose();
  qTransform = tAdjoint * qTransform;

  setCoef(2, 0, 0,     qTransform[0][0]);
  setCoef(1, 1, 0, 2.0*qTransform[0][1]);
  setCoef(1, 0, 1, 2.0*qTransform[0][2]);
  setCoef(1, 0, 0, 2.0*qTransform[0][3]);
  setCoef(0, 2, 0,     qTransform[1][1]);
  setCoef(0, 1, 1, 2.0*qTransform[1][2]);
  setCoef(0, 1, 0, 2.0*qTransform[1][3]);
  setCoef(0, 0, 2,     qTransform[2][2]);
  setCoef(0, 0, 1, 2.0*qTransform[2][3]);
  setCoef(0, 0, 0,     qTransform[3][3]);
}

gmVector3 Quadric::Centroid()
{
  return gmVector3(0,0,0);
}

