#ifndef __CONE_H
#define __CONE_H

#include "Quadric.h"

class Cone : public Quadric
{
  public:
    Cone() : Quadric(1.0, 1.0, -1.0, 0.0) { }          // apex at origin
    Cone(double r) : Quadric(1.0, 1.0, -r*r, 0.0) { }  // radius one from apex
    Cone(const gmVector3& d, double r);                       // direction of cone axis
    Cone(const gmVector3 & x, gmVector3 d, double r);          // position of apex

    MAKE_NAME();

    virtual const char ** getPixmapXPM(const int& size) const;
};

#endif

