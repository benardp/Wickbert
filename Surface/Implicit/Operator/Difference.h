/** 
 * @file Difference.h 
 * Definition of the Difference implicit class
 * @author John C. Hart
 * @date 4 Oct. 2001
 */

#ifndef DIFFERENCE_H
#define DIFFERENCE_H

#include "BinaryOp.h"

/** Difference is a class of Implicit that represents the sum of two
 * implicits.
 */
class Difference : public BinaryOp
{
  public:
    /** Constructor
     * Defaults to binary op constructor.
     */
    Difference() { m_f = m_g = NULL; };

#ifndef INTERVAL_EVAL_ONLY
    virtual double proc(const gmVector3 & x);
    virtual gmVector3 grad(const gmVector3 & x);
    virtual gmMatrix3 hess(const gmVector3 & x);
#endif

    virtual Intervald proc(const Box<double>& );
    virtual Box3d grad(const Box<double>& );
    virtual IMatrix3d hess(const Box<double>& );

    MAKE_NAME();
};

#endif

