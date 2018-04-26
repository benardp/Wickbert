/**
* Declaration of a binary operation.
* @file BinaryOp.h
* @date Fall 200
* @author John Hart, Ed Bachta
* @remarks Based on UnaryOp.h by Ed Bachta and BinaryOp by Bill Nagel
*/

#ifndef __ASL_BINARYOP_H_
#define __ASL_BINARYOP_H_

#include "Surface/Implicit/Implicit.h"

class BinaryOp : public Implicit
{
  public:
    Implicit *m_f; ///< First operand.
    Implicit *m_g; ///< Second operand.

    BinaryOp(Implicit *f = NULL, Implicit *g = NULL); ///< Explicit constructor.

    virtual unsigned int qlen();
    virtual void _setq(double *q);
    virtual void getq(double *q);
    virtual void procq(const gmVector3 & x, double *q);

    virtual void getqname(char **qn);
    virtual bool setChild(int index, Implicit* child);
    virtual Implicit* getChild(int index);

    virtual int maxChildren() { return 2; };
    virtual int numChildren() ;
};

#endif

