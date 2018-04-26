#include "Torus.h"

/** Add Torus to the ImplicitFactory
 */
REGISTER_IMPLICIT(Torus,"Algebraic:Torus");

Torus::Torus() : Algebraic(4)
{
  setCoef(4,0,0,1.0);
  setCoef(0,4,0,1.0);
  setCoef(0,0,4,1.0);
  setCoef(2,2,0,2.0);
  setCoef(0,2,2,2.0);
  setCoef(2,0,2,2.0);
  setCoef(2,0,0,-10.0);
  setCoef(0,2,0,-10.0);
  setCoef(0,0,2,6.0);
  setCoef(0,0,0,9.0);
}

Torus::Torus(double R, double r) : Algebraic(4)
{
  double R2 = R*R;
  double r2 = r*r;
  double R4 = R2*R2;
  double r4 = r2*r2;

  setCoef(4,0,0,1.0);
  setCoef(0,4,0,1.0);
  setCoef(0,0,4,1.0);
  setCoef(2,2,0,2.0);
  setCoef(0,2,2,2.0);
  setCoef(2,0,2,2.0);
  setCoef(2,0,0,-2.0*(R2 + r2));
  setCoef(0,2,0,-2.0*(R2 + r2));
  setCoef(0,0,2,2.0*(R2 - r2));
  setCoef(0,0,0,R4 + r4 - 2*R2*r2);
}

const char ** Torus::getPixmapXPM(const int& size) const
{
  if (size <= 16)
    return (const char **)torus_pixmap16;
  else if (size <= 32)
    return (const char **)torus_pixmap32;
  else
    return (const char **)torus_pixmap48;
}

