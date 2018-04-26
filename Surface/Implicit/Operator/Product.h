/** 
 * @file Product.h 
 * Definition of the Product implicit operator
 * @author John C. Hart
 * @date 6 Dec. 2005
 */

#ifndef PRODUCT_H
#define PRODUCT_H

#include "BinaryOp.h"

/** Product is a class of Implicit that represents the product of two implicits.
 */
class Product : public BinaryOp
{
public:
    Product(Implicit *f = NULL, Implicit *g = NULL) { m_f = f; m_g = g; }

    virtual double proc(const gmVector3 & x);
    virtual gmVector3 grad(const gmVector3 & x);
    virtual gmMatrix3 hess(const gmVector3 & x);

    virtual Intervald proc(const Box<double>&  b);
    virtual Box3d grad(const Box<double>&  b);
    virtual IMatrix3d hess(const Box<double>&  b);

    MAKE_NAME();
};

#endif 

