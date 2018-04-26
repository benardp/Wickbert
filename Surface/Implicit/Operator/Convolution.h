/** @file Convolution.h
 * Sums implicits placed at particles
 * @author: John C. Hart
 * @date: 10 Jan. 2006
 */

#ifndef CONVOLUTION_H
#define CONVOLUTION_H

#include "UnaryOp.h"
#include "Particles/Attributes/ParticlePosition.h"

/** Returns the sum of copies of the implicit operand at each of the particle positions
 */
class Convolution : public UnaryOp
{
public:
	ParticlePosition *positions;
	  
    /// Constructor.
    Convolution(void);

	virtual double proc(const gmVector3 & x) {
		double sum = 0.0;
		if (m_f && positions)
			for (unsigned int i = 0; i < positions->x.size(); i++)
				sum += m_f->proc(x - positions->x[i]);
		return sum;
	}

	virtual gmVector3 grad(const gmVector3 & x) {
		gmVector3 sum;
		if (m_f && positions)
			for (unsigned int i = 0; i < positions->x.size(); i++)
				sum += m_f->grad(x - positions->x[i]);
		return sum;
	}

	virtual gmMatrix3 hess(const gmVector3 & x) {
		gmMatrix3 sum;
		if (m_f && positions)
			for (unsigned int i = 0; i < positions->x.size(); i++)
				sum += m_f->hess(x - positions->x[i]);
		return sum;
	}
	
	/// Functions for dealing with the parameter vector.
	virtual unsigned int qlen() { return (m_f ? m_f->qlen() : 0) + (positions ? 3*positions->x.size() : 0); }

	virtual void getq(double *q) {
		int n = m_f ? m_f->qlen() : 0;
		if (m_f) m_f->getq(q); 
		if (positions)
			for (unsigned int i = 0; i < positions->x.size(); i++)
				positions->x[i].copyTo(&q[n + 3*i]);
	}

	virtual void getqname(char **qn) {
		int n = m_f ? m_f->qlen() : 0;
		if (m_f) m_f->getqname(qn);
		if (positions) {
			for (unsigned int i = 0; i < positions->x.size(); i++) {
				qn[n + 3*i] = "x";
				qn[n + 3*i + 1] = "y";
				qn[n + 3*i + 2] = "z";
			}
		}
	}

	virtual void _setq(double *q) {
		int n = m_f ? m_f->qlen() : 0;
		if (m_f) m_f->setq(q);
		if (positions)
			for (unsigned int i = 0; i < positions->x.size(); i++)
				positions->x[i].assign(q[n + 3*i],q[n + 3*i + 1],q[n + 3*i + 2]);
	}
    
    MAKE_NAME();
};

#endif 

