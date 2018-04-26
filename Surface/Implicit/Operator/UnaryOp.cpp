/**
 * Implementation of the unary operation.
 * @file UnaryOp.cpp
 * @date July 13, 2001
 * @author Ed Bachta
 */

#include "UnaryOp.h"

UnaryOp::UnaryOp(Implicit *f)
{
	new SurfImpRefParam(this,&m_f);
}


/**
 * Sets one of the operands of this operation.
 * @param   index the operand to set (either 0, or 1)
 * @param   child What to set the operand to.
 * @returns False if there are not exactly 2 children in the vector.
 */
bool UnaryOp::setChild(int index, Implicit* child)
{
  bool retval = false;

  if (index == 0)
    {
      m_f = child;
      retval = true;
    }
  
  return retval;
}

/**
 * Gets one of the operands of this operation.
 * @param   index the operand to get (either 0, or 1)
 * @returns child[index].
 */
Implicit* UnaryOp::getChild(int index)
{
  if (index == 0)
    return m_f;

  return NULL;
}

int UnaryOp::numChildren()
{
  int numchild = 0;
  if (m_f != NULL)
    numchild++;

  return numchild;
}

void UnaryOp::getqname(char **qn)
{
  int myqlen;
  std::string name;

  if (!m_f) 
    return;

  myqlen = qlen() - m_f->qlen();

  m_f->getqname(&qn[myqlen]);

#ifdef USENAMES
  // Insert the subobject names in parameter names
  for (int i = myqlen; i < myqlen + m_f->qlen(); i++)
    {
      name = m_f->getObjectName() + ":" + qn[i];
      qn[i] = (char *)malloc(name.length()*sizeof(char) + 1);
      strcpy(qn[i],name.c_str());
  }
#endif
}

