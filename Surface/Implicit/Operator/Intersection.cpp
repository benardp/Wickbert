#include "Intersection.h"

#if 0
REGISTER_IMPLICIT(Intersection,"BinaryOp:Intersection");

/**
 * Called by the various constructors to allow for a single location for
 * the init of an Intersection object.
 * @param f    One function to be intersected.
 * @param g    The other function to be intersected.
 * @param cont  The continuity to be provided.
 */
void Intersection::init(Implicit* f, Implicit* g, double cont)
{
  m_f = f;
  m_g = g;
  m_cont = cont;
  m_sign = 1;
}

/**
 * Default constructor.
 */
Intersection::Intersection() 
{
  init(NULL,NULL,1.0);
}

/**
 * Constructor for an Intersection.
 * @param f    One function to be intersected.
 * @param g    The other function to be intersected.
 */
Intersection::Intersection(Implicit *f, Implicit *g) 
{ 
  init(f,g,1.0);
}

/**
 * Constructor for an Intersection.
 * @param f    One function to be intersected.
 * @param g    The other function to be intersected.
 * @param cont  The continuity to be provided.
 */
Intersection::Intersection(Implicit *f, Implicit *g, double cont) 
{ 
  init(f,g,cont);
}
#endif

REGISTER_IMPLICIT(Intersection, "BinaryOp:Intersection");

