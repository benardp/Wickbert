/** \file Quadric.h Specification of a quadric surface class
 * \author Jeff Decker
 * \date Fall 2000
 */

#ifndef __QUADRIC_H
#define __QUADRIC_H

#include "Algebraic.h"

/** 
 * A quadric implicit surface class.
 * The Quadric class is a special subset of the Algebraic class with added
 * functionality that can take advantage of the homogeneous respresentation
 * of quadrics.
 */
class Quadric : public Algebraic
{
  public:
    Quadric() : Algebraic(2) { } 
    Quadric(double a, double e, double h, double j); 
    Quadric(double a, double b, double c, double d, double e, 
            double f, double g, double h, double i, double j); 
    Quadric(double q[10]); 

    gmMatrix4 Q(); 
    void Transform(gmMatrix4 t); 
    gmVector3 Centroid();

    MAKE_NAME();
};

#endif

