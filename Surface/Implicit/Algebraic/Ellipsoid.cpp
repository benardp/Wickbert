#include "Ellipsoid.h"

// Add Ellipsoid to the ImplicitFactory
REGISTER_IMPLICIT(Ellipsoid,"Algebraic:Quadric:Ellipsoid");

Ellipsoid::Ellipsoid(const gmVector3 & x, double a, double b, double c) : 
Quadric(1.0/(a*a), 1.0/(b*b), 1.0/(c*c), -1.0)
{
  gmMatrix4 transMatrix = gmMatrix4::translate(x[0], x[1], x[2]);

  this->Transform(transMatrix);
}

const char ** Ellipsoid::getPixmapXPM(const int& size) const
{
  if (size <= 16)
    return (const char **)ellipsoid_pixmap16;
  else if (size <= 32)
    return (const char **)ellipsoid_pixmap32;
  else
    return (const char **)ellipsoid_pixmap48;
}

