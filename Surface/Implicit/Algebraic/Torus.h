#ifndef __TORUS_H
#define __TORUS_H

#include "Algebraic.h"

class Torus : public Algebraic
{
  public:
    Torus();
    Torus(double R, double r);  // torus with major radius R and minor radius r
                                // surrounding z-axis 

    MAKE_NAME();

    virtual const char ** getPixmapXPM(const int& size) const;
};

#endif

