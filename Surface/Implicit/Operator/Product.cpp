/**
 * @file Product.cpp
 * Implementation of product operation
 */
#include "Product.h"

REGISTER_IMPLICIT(Product,"BinaryOp:Product");

double Product::proc(const gmVector3 & x) 
{
	return (m_f && m_g) ? m_f->proc(x)*m_g->proc(x) : 0.0;
}

/** Implements product rule for differentiation
 */
gmVector3 Product::grad(const gmVector3 & x) 
{
	return (m_f && m_g) ? m_f->proc(x)*m_g->grad(x) + m_f->grad(x)*m_g->proc(x) : gmVector3();
}

/** (f*g)'' = (f'*g + f*g')' = f''*g + outer(f',g') + outer(g',f') + f*g''
 */
gmMatrix3 Product::hess(const gmVector3 & x) 
{
	return (m_f && m_g) ? m_f->hess(x)*m_g->proc(x) + outer(m_f->grad(x),m_g->grad(x)) +
						  outer(m_g->grad(x),m_f->grad(x)) + m_f->proc(x)*m_g->hess(x)
						: gmMatrix3();
}

Intervald Product::proc(const Box<double>&  b) 
{ 
	return (m_f && m_g) ? m_f->proc(b)*m_g->proc(b) : Intervald(0.0);
}

Box3d Product::grad(const Box<double>&  b) 
{ 
	//stupid gcc 3.3 needs this _temp --ms
	Box3d _temp;
	if(m_f && m_g)
	 _temp =  m_f->proc(b)*m_g->grad(b) + m_f->grad(b)*m_g->proc(b);
	
	return _temp;
}

IMatrix3d Product::hess(const Box<double>&  b) 
{ 
	//stupid gcc 3.3 needs this _temp --ms
	IMatrix3d _temp;
	if(m_f && m_g)
	_temp =  m_f->hess(b)*m_g->proc(b) + outer(m_f->grad(b),m_g->grad(b)) + outer(m_g->grad(b),m_f->grad(b)) + m_f->proc(b)*m_g->hess(b);
	
	return _temp; 
}

