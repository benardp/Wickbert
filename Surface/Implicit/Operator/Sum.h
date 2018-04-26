/** 
 * @file Sum.h 
 * Definition of the Sum implicit class
 * @author John C. Hart
 * @date 25 Sep. 2001
 * 
 * Implements the Sum implicit class.
 */

#ifndef SUM_H
#define SUM_H

#include "BinaryOp.h"

/** Sum is a class of Implicit that represents the sum of two implicits.
 */
class Sum : public BinaryOp
{
  public:
    /** Constructor
     * Defaults to binary op constructor.
     */
    Sum() { m_f = m_g = NULL; }

#ifndef INTERVAL_EVAL_ONLY
    virtual double proc(const gmVector3 & x);
    virtual gmVector3 grad(const gmVector3 & x);
    virtual gmMatrix3 hess(const gmVector3 & x);
#endif

    virtual Intervald proc(const Box<double>&  b);
    virtual Box3d grad(const Box<double>&  b);
    virtual IMatrix3d hess(const Box<double>&  b);

    MAKE_NAME();
};

#endif 

