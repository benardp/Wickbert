/**
 * Declaration of the unary operation.
 * @file UnaryOp.h
 * @date Fall 2000
 * @author Ed Bachta
 */

#ifndef __ASL_UNARYOP_H__
#define __ASL_UNARYOP_H__

#include "Surface/Implicit/Implicit.h"

class UnaryOp : public Implicit
{
public:
  Implicit *m_f; ///< Surface on which to operate

  /** Default constructor.
   * Sets m_f to NULL.
   */
  UnaryOp(Implicit *f = NULL);

  /** qlen Automatically returns total qlen of children.
   * Does not need to be overridden if op has no parameters.
   * Can be called from subclasses to handle children.
   */
  virtual unsigned int qlen() { return m_f ? m_f->qlen() : 0; }

  /** _setq Automatically sets q of children.
   * Does not need to be overridden if op has no parameters.
   * Can be called from subclasses to handle children.
   */
  virtual void _setq(double *q) { if (m_f) m_f->_setq(q); }

  /** getq Automatically gets q of children.
   * Does not need to be overridden if op has no parameters.
   * Can be called from subclasses to handle children.
   */
  virtual void getq(double *q) { if (m_f) m_f->getq(q); }

#ifndef INTERVAL_EVAL_ONLY
  virtual double proc(const gmVector3 & x)
    { return m_f ? m_f->proc(x) : 0.0; }
  virtual gmVector3 grad(const gmVector3 & x) 
    { return m_f ? m_f->grad(x) : gmVector3(); }
  virtual gmMatrix3 hess(const gmVector3 & x)
    { return m_f ? m_f->hess(x) : gmMatrix3(); }
#endif

  virtual Intervald proc(const Box<double>& x)
    { return m_f ? m_f->proc(x) : Intervald(0.0); }
  virtual Box3d grad(const Box<double>& x)
    { return m_f ? m_f->grad(x) : Box3d(0.0); }
  virtual IMatrix3d hess(const Box<double>& x)
    { return m_f ? m_f->hess(x) : IMatrix3d(0.0); }

  /** procq Automatically computes df/dq of children.
   * Does not need to be overridden if op has no parameters.
   * Can be called from subclasses to handle children.
   */
  virtual void procq(const gmVector3 & x, double *q) 
    //{ if (m_f) m_f->procq(x, &q[qlen()]); }
    { if (m_f) m_f->procq(x, q); }

  /** Automatically fills qn with operand parameter names.
   * UnaryOp's with no parameters need not redefine getqname().
   * UnaryOp's with parameters should set the names of only their parameters
   * and then call UnaryOp::getqname(qn) to let it set its operand's
   * parameters.
   */
  virtual void getqname(char **qn);

  /// Sets one of the operands of this operation.
  virtual bool setChild(int index, Implicit* child);

  /// Get the specified child
  virtual Implicit* getChild(int index);

  virtual int maxChildren() { return 1; };

  virtual int numChildren() ;
};

#endif

