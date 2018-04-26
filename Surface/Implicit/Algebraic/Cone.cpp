#include "Cone.h"

/** Registers Cone in the Implicit factory.
 */
REGISTER_IMPLICIT(Cone,"Algebraic:Quadric:Cone");

Cone::Cone(const gmVector3& d, double r) : Quadric(1.0, 1.0, -r*r, 0.0)
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
  Transform(rotMatrix);
}

Cone::Cone(const gmVector3 & x, gmVector3 d, double r) : Quadric(1.0, 1.0, -r*r, 0.0)
{
  gmVector3 zAxis, rotAxis;
  gmMatrix4 rotMatrix, transMatrix;
  double rotCosine;
  double rotAngle;

  zAxis.assign(0.0, 0.0, 1.0);
  rotAxis = cross(d, zAxis);
  rotCosine = dot(d, zAxis); // d is already normalized
  rotAngle = acos(rotCosine);
  rotMatrix = gmMatrix4::rotate(rotAngle, rotAxis);
  Transform(rotMatrix);

  transMatrix = gmMatrix4::translate(x[0],x[1],x[2]);
  Transform(transMatrix);
}

const char ** Cone::getPixmapXPM(const int& size) const
{
  if (size <= 16)
    return (const char **)cone_pixmap16;
  else if (size <= 32)
    return (const char **)cone_pixmap32;
  else
    return (const char **)cone_pixmap48;
}

