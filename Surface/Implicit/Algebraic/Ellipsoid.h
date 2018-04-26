#ifndef __ELLIPSOID_H
#define __ELLIPSOID_H

#include "Quadric.h"

class Ellipsoid : public Quadric
{
  public:
    Ellipsoid() : Quadric(1.0, 1.0, 1.0, -1.0) { }
    Ellipsoid(double a, double b, double c) :
      Quadric(1.0/(a*a), 1.0/(b*b), 1.0/(c*c), -1.0) { }
    Ellipsoid(const gmVector3 & x, double a, double b, double c);  // position

    MAKE_NAME();

    virtual const char ** getPixmapXPM(const int& size) const;
};

#endif

