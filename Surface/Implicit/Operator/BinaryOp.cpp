/**
 * Implementation of the binary operation.
 * @file UnaryOp.cpp
 * @date July 13, 2001
 * @author Ed Bachta
 */

#include "BinaryOp.h"

/**
 * Explicit constructor.
 * @param f First operand.
 * @param g Second operand.
 */
BinaryOp::BinaryOp(Implicit *f, Implicit *g)
{
  m_f = f;
  m_g = g;

  new SurfImpRefParam(this,&m_f);
  new SurfImpRefParam(this,&m_g);
} 

/** 
 * qlen Automatically returns total qlen of children.
 * Does not need to be overridden if op has no parameters.
 * Can be called from subclasses to handle children.
 */
unsigned int BinaryOp::qlen() 
{ 
  return ((m_f ? m_f->qlen() : 0) + 
          (m_g ? m_g->qlen() : 0)); 
}

/** 
 * setq Automatically sets q of children.
 * Does not need to be overridden if op has no parameters.
 * Can be called from subclasses to handle children.
 */
void BinaryOp::_setq(double *q) 
{
  if (m_f) m_f->_setq(q);
  if (m_g) m_g->_setq(&q[(m_f ? m_f->qlen() : 0)]);
}

/**
 * getq Automatically gets q of children.
 * Does not need to be overridden if op has no parameters.
 * Can be called from subclasses to handle children.
 */
void BinaryOp::getq(double *q) 
{
  if (m_f) m_f->getq(q);
  if (m_g) m_g->getq(&q[(m_f ? m_f->qlen() : 0)]);
}

/** 
 * procq Automatically computes df/dq of children.
 * Does not need to be overridden if op has no parameters.
 * Can be called from subclasses to handle children.
 */
void BinaryOp::procq(const gmVector3 & x, double *q) 
{
  if (m_f) m_f->procq(x, q);
  if (m_g) m_g->procq(x, &q[(m_f ? m_f->qlen() : 0)]);
}

/**
 * Sets one of the operands of this operation.
 * @param   index the operand to set (either 0, or 1)
 * @param   child What to set the operand to.
 * @returns False if there are not exactly 2 children in the vector.
 */
bool BinaryOp::setChild(int index, Implicit* child)
{
  bool retval = false;

  if (index == 0)
    {
      m_f = child;
      retval = true;
    }
  else if (index == 1)
    {
      m_g = child;
      retval = true;
    }
  
  return retval;
}

/**
* Gets one of the operands of this operation.
* @param   index the operand to get (either 0, or 1)
* @returns child[index].
*/
Implicit* BinaryOp::getChild(int index)
{
  if (index == 0)
    return m_f;
  else if (index == 1)
    return m_g;
  
  return NULL;
}

int BinaryOp::numChildren()
{
  int numchild = 0;
  if (m_f != NULL)
    numchild++;
  if (m_g != NULL)
    numchild++;

  return numchild;
}

/** 
 * Automatically fills qn with operands parameter names.
 * BinaryOp's with no parameters need not redefine getqname().
 * BinaryOp's with parameters should set the names of only their parameters
 * and then call BinaryOp::getqname(qn) to let it set its operands'
 * parameters.
 */
void BinaryOp::getqname(char **qn)
{
  std::string name;
  int fqlen = m_f ? m_f->qlen() : 0;
  int gqlen = m_g ? m_g->qlen() : 0;
  int myqlen = qlen() - fqlen - gqlen;
  
  if (m_f) m_f->getqname(&qn[myqlen]);
  if (m_g) m_g->getqname(&qn[myqlen + fqlen]);
  
#ifdef USENAMES
  // Insert the subobject names in parameter names
  for (i = myqlen; i < myqlen + fqlen; i++)
    {
      name = m_f->getObjectName() + ":" + qn[i];
      qn[i] = (char *)malloc(name.length()*sizeof(char) + 1);
      strcpy(qn[i],name.c_str());
    }
  
  for (i = myqlen + fqlen; i < myqlen + fqlen + gqlen; i++)
    {
      name = m_g->getObjectName() + ":" + qn[i];
      qn[i] = (char *)malloc(name.length()*sizeof(char) + 1);
      strcpy(qn[i],name.c_str());
    }
#endif
}

