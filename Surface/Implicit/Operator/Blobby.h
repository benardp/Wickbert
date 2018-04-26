/** @file Convolution.h
 * Sums implicits placed at particles
 * @author: John C. Hart
 * @date: 10 Jan. 2006
 */

#ifndef BLOBBY_H
#define BLOBBY_H

#include "Particles/Attributes/ParticlePosition.h"
#include "Particles/Attributes/ParticleScalar.h"
#include "../Blob/Blinn.h"
#include "Implicit/Algebraic/Quadric.h"

/** Returns the sum of copies of the implicit operand at each of the particle positions
 */
class Blobby : public Implicit
{
public:
	ParticlePosition *positions;
	ParticleScalar *radii;
	Quadric *basis;
	Blinn *blinn;

	double blobbiness;

    /// Constructor.
    Blobby(void);
   
	virtual double proc(const gmVector3 & x) {
		blinn->m_b = blobbiness;
		double sum = 0.0;
		if (positions && radii) {
			for (unsigned int i = 0; i < positions->x.size(); i++) {
				blinn->m_r = radii->getScalar(i);
				sum += blinn->proc(x - positions->x[i]);
			}
		}
		return 0.5 - sum;
	}

	virtual gmVector3 grad(const gmVector3 & x) {
		blinn->m_b = blobbiness;
		gmVector3 sum;
		if (positions && radii) {
			for (unsigned int i = 0; i < positions->x.size(); i++) {
				blinn->m_r = radii->getScalar(i);
				sum += blinn->grad(x - positions->x[i]);
			}
		}
		return -sum;
	}

	virtual gmMatrix3 hess(const gmVector3 & x) {
		blinn->m_b = blobbiness;
		gmMatrix3 sum;
		if (positions && radii) {
			for (unsigned int i = 0; i < positions->x.size(); i++) {
				blinn->m_r = radii->getScalar(i);
				sum += blinn->hess(x - positions->x[i]);
			}
		}
		return -sum;
	}
	
	// Functions for dealing with the parameter vector.

	virtual unsigned int qlen() { return (positions && radii)? 4*positions->x.size() : 0; }

	virtual void getq(double *q) {
		if (positions && radii) {
			for (unsigned int i = 0; i < positions->x.size(); i++) {
				positions->x[i].copyTo(&q[4*i]);
				q[4*i + 3] = radii->getScalar(i);
			}
		}
	}

	virtual void getqname(char **qn) {
		if (positions && radii) {
			for (unsigned int i = 0; i < positions->x.size(); i++) {
				qn[4*i] = "x";
				qn[4*i + 1] = "y";
				qn[4*i + 2] = "z";
				qn[4*i + 3] = "r";
			}
		}
	}

	virtual void _setq(double *q) {
		if (positions && radii) {
			for (unsigned int i = 0; i < positions->x.size(); i++) {
				positions->x[i].assign(q[4*i],q[4*i + 1],q[4*i + 2]);
				radii->setScalar(i,q[4*i + 3]);
			}
		}
	}
    
    MAKE_NAME();
};

#endif 

