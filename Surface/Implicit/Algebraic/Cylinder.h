#ifndef __CYLINDER_H
#define __CYLINDER_H

#include "Quadric.h"

class Cylinder : public Quadric
{
  public:
    Cylinder() : Quadric(1.0, 1.0, 0.0, -1.0) { }
    Cylinder(double r) : Quadric(1.0, 1.0, 0.0, -r*r) { }
    Cylinder(const gmVector3& d, double r);
    Cylinder(const gmVector3 & x, gmVector3 d, double r);

    MAKE_NAME();

    virtual const char ** getPixmapXPM(const int& size) const;
}; 

#endif

