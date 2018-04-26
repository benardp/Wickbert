#include "Cylinder.h"

// Add Cylinder to the ImplicitFactory
REGISTER_IMPLICIT(Cylinder,"Algebraic:Quadric:Cylinder");

Cylinder::Cylinder(const gmVector3& d, double r) : Quadric(1.0, 1.0, 0.0, -r*r)
{
  gmVector3 zAxis, rotAxis;
  gmMatrix4 rotMatrix;
  double rotCosine;
  double rotAngle;

  zAxis.assign(0.0, 0.0, 1.0);
  rotAxis = cross(d, zAxis);
  rotCosine = dot(d, zAxis);  // d is already normalized
  rotAngle = acos(rotCosine);
  rotMatrix = gmMatrix4::rotate(rotAngle, rotAxis);
  this->Transform(rotMatrix);
}

Cylinder::Cylinder(const gmVector3 & x, gmVector3 d, double r) : 
Quadric(1.0, 1.0, 0.0, -r*r)
{
  gmVector3 zAxis, rotAxis;
  gmMatrix4 rotMatrix, transMatrix;

/* SOMETHING WRONG WITH THIS!
  double rotCosine;
  double rotAngle;
  zAxis.assign(0.0, 0.0, 1.0);
  rotAxis = cross(d, zAxis);
  rotCosine = dot(d, zAxis);  // d is already normalized
  rotAngle = acos(rotCosine);
  rotMatrix = gmMatrix4::rotate(rotAngle, rotAxis);
  this->Transform(rotMatrix);
*/
  transMatrix = gmMatrix4::translate(x[0], x[1], x[2]);
  this->Transform(transMatrix);
}

const char ** Cylinder::getPixmapXPM(const int& size) const
{
  if (size <= 16)
    return (const char **)cylinder_pixmap16;
  else if (size <= 32)
    return (const char **)cylinder_pixmap32;
  else
    return (const char **)cylinder_pixmap48;
}

