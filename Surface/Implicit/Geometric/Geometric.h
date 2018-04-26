/**
 * @file Geometric.h
 */
#ifndef __GEOMETRIC_H
#define __GEOMETRIC_H

#include "Surface/Implicit/Implicit.h"

/** Class of implicit surface primitives based on distance.
 *  The geometric class holds all geometric methods dealing with implicit
 *  surfaces based on distance.
 */
class Geometric : public Implicit
{
public:
  Geometric() : Implicit() { }

  /** Gradient of a function returning distance is everywhere magnitude one.
   */
  virtual double lipschitz(const Box<double>& domain) { return 1.0; }
};

#endif 

