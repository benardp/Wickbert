/**
* Implementation of the implicit list.
* @file ImpList.cpp
* @date July 18, 2001
* @author Ed Bachta
*/

#include "ImpList.h"

REGISTER_IMPLICIT(ImpList,"ImpList");

/**
 * Default constructor.
 */
ImpList::ImpList() 
{
}

/**
 * Explicit constructor.
 */
ImpList::ImpList(std::list<Implicit*> ilist) : std::list<Implicit*>(ilist) 
{
}

#ifndef INTERVAL_EVAL_ONLY
/**
 * Evaluate the sum of f(x).
 * @param   x Point at which to evaluate.
 * @returns Sum of f(x) over implicits in the list.
 */
double ImpList::proc(const gmVector3 & x)
{
  double result = 0.0;
  
  for (i_it = this->begin(); i_it != this->end(); i_it++) 
    result += (*i_it)->proc(x);
  
  return result;
} 

/**
 * Evaluate the sum of grad(x).
 * @param   x Point at which to evaluate.
 * @returns Sum of grad(x) over implicits in the list.
 */
gmVector3 ImpList::grad(const gmVector3 & x) 
{
  gmVector3 result;
  
  for (i_it = this->begin(); i_it != this->end(); i_it++) 
    result += (*i_it)->grad(x);
  
  return result;
}

gmMatrix3 ImpList::hess(const gmVector3 & x) 
{
  gmMatrix3 result;
  
  for (i_it = this->begin(); i_it != this->end(); i_it++) 
    result = result + (*i_it)->hess(x);
  
  return result;
}
#endif

Intervald ImpList::proc(const Box<double>& x)
{
  Intervald result = Intervald(0.0);
  
  for (i_it = this->begin(); i_it != this->end(); i_it++) 
    result += (*i_it)->proc(x);
  
  return result;
} 

Box3d ImpList::grad(const Box<double>& x) 
{
  Box3d result = Box3d(0.0);
  
  for (i_it = this->begin(); i_it != this->end(); i_it++) 
    result = result + (*i_it)->grad(x);
  
  return result;
}

IMatrix3d ImpList::hess(const Box<double>& x) 
{
  IMatrix3d result;
  
  for (i_it = this->begin(); i_it != this->end(); i_it++) 
    result = result + (*i_it)->hess(x);
  
  return result;
}

/**
 * Return number of parameters.
 */
unsigned int ImpList::qlen() 
{
  int result = 0;
  
  for (i_it = this->begin(); i_it != this->end(); i_it++) 
    result += (*i_it)->qlen();
  
  return result;
} // end qlen

/**
 * Retrieve parameters.
 * @param q Array of parameters.
 */
void ImpList::getq(double* q) 
{
  int pos = 0;
  for (i_it = this->begin(); i_it != this->end(); i_it++) 
    {
      (*i_it)->getq(&q[pos]);
      pos += (*i_it)->qlen();
    }
} 

/**
 * Assign parameters.
 * @param q Array of parameters.
 */
void ImpList::_setq(double* q) 
{
  int pos = 0;
  for (i_it = this->begin(); i_it != this->end(); i_it++) 
    {
      (*i_it)->_setq(&q[pos]);
      pos += (*i_it)->qlen();
    }
}

/**
 * Retrieve parameter names.
 * @param n Array of parameter names.
 */
void ImpList::getqname(char** n)
{  
  int pos = 0;
  for (i_it = this->begin(); i_it != this->end(); i_it++)
    {    
      (*i_it)->getqname(&n[pos]);
      pos += (*i_it)->qlen();    
    }
}

/**
 * Assign elements.
 * @param index The index to place the child.
 * @param child The child to put into the ImpList.
 */
bool ImpList::setChild(int index, Implicit* child)
{
  if (index < (int)this->size())
    {  // replace an element
      i_it = this->begin();
      std::advance(i_it, index);
      (*i_it) = child;
      return true;
    }
  else if (index == this->size())
    {  // insert an element
      push_back(child);
      return true;
    }
  else  // index out of range
    return false;
}

/**
 *  Retrieve elements.
 * @index The index of the child to be retrieved.
 */
Implicit* ImpList::getChild(int index)
{
  i_it = this->begin();
  std::advance(i_it, index);
  return *i_it;
}

