/**
 * Represents the union of two implicit surfaces.
 *
 * @author Steve Zelinka
 */

#ifndef _UNION_H
#define _UNION_H

#include "RFunction.h"

/** Union creates an Rfunction with default C1 continuity.
 */
class Union : public RFunction
{
  public:
    /**
    * Constructor for a Union.
    *
    * @param f    One function to be union'd.
    * @param g    The other function to be union'd.
    * @param cont  The degree of continuity to be provided.
    */
	 Union(Implicit *f, Implicit *g, int cont = 1) : RFunction(cont, -1.0)
    { 
      m_f = f; 
      m_g = g; 
      //m_cont = cont;
      //m_sign = -1;
    }

	Union() :RFunction(1, -1) {
      m_f = m_g = NULL;
      //m_cont = 1;
      //m_sign = -1;
    }

    MAKE_NAME();
};

#endif

