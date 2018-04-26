/**
* Implementation of SurfaceDeformation.
* @file SurfaceDeformation.cpp
* @author John Hart, Ed Bachta
*/

#include "SurfaceDeformation.h"

#include "Surface/Implicit/Implicit.h"
#include "Surface/Implicit/Variational/RBF.h"

#include "ImplicitInterrogator.h"
#include "ParticlePosition.h"
#include "ParticleVelocity.h"
#include "ParticleScalar.h"
#include "ParticleNormal.h"

// TNT includes
#include "tnt/tnt.h"
#include "tnt/vec.h"
#include "tnt/cholesky.h"
#include "tnt/cmat.h"
#include "tnt/transv.h"
#include "tnt/trisolve.h"

REGISTER_PARTICLESTUFF(SurfaceDeformation,"Behavior:SurfaceDeformation");

/**
* Creates a surface deformation behavior for Particles p.
* @param ps The owning particle system.
*
* Sets phi (the feedback constant) to 15.
* 
*/
SurfaceDeformation::SurfaceDeformation(Particles *ps)
: ParticleBehavior(ps, std::string("SurfaceDeformation"))
{
	new PSParamDouble(this,&phi,15.0,"k","Feedback",
		"Strength of penalty force to keep surface through particle.");

	new Attached<ImplicitInterrogator>(this,&imp_int);
	new Attached<ParticleNormal>(this,&p_orient);
	new Attached<ParticlePosition>(this,&position);
	new Attached<ParticleVelocity>(this,&velocity);
	new Attached<ParticleScalar>(this,&scalar);
}

void SurfaceDeformation::attachAttributes()
{
	ParticleBehavior::attachAttributes();

	if (imp_int)
	{
		// set the flexible length to the size of imp_int
		Implicit* imp=imp_int->getImplicit();
		if (imp)
		{
			flexible.resize(imp->qlen());
			for(unsigned int i=0;i<imp->qlen();i++)
				flexible[i]=true;
		}
	}
}

void SurfaceDeformation::loadParticles()
{
	unsigned int i;

	if (!imp_int)
		return;
	RBF *rbf = dynamic_cast<RBF *>(imp_int->getImplicit());
	
	if(!rbf)
		return;

	else
	{
		for (i = 0; i < rbf->constraints.size(); i++)
		{
			ps->addParticle();
			position->setPosition(i,rbf->constraints[i].c);
			h[i] = rbf->constraints[i].h;
			//d_orient->n[i] = rbf->constraints[i].n;
		}
	}
}

void SurfaceDeformation::particleAdded()
{
	position->changed = true;
	h.push_back(0.0);
}

void SurfaceDeformation::particleRemoved(unsigned int i)
{
	position->changed = true;
	h[i] = h.back();
	h.pop_back();
}

void SurfaceDeformation::applyConstraint()
{
	unsigned int i;

	if (!imp_int) return;

	if ( ps->findAttribute(position) && position->changed ) {
		Particles *tp = ps;   //ps->particleSystem->findParticles(srcparticles);
		if ( tp )
			imp_int->getImplicit()->interpolate(tp, (float) phi);
	}

	for (i=0; i<ps->size(); i++)
	{
		gmVector3 gradient = imp_int->grad(i);//getImplicit()->grad(position->getPosition(i));
		
		p_orient->setNormal(i, gradient);
		p_orient->setGradMag2(i, gradient.lengthSquared());
	}
}

/// Change the implicit to pass through particles
// [tk 5.19.05] function should be moved to Implicit
// This default behavior should be a virtual function in implicit, which should be overridden in RBF and other derived classes that need special functionality
// I think this is how the extra ParticleScalar object gets loaded into the default particles system - it's there because the SurfaceDeformation might need it to connect
// to an RBF.  This should be unnecessary for most uses of SurfaceDeformation, and it should be removable once we move the functionality from interpolate(...) to its proper
// location.
#if 0
void SurfaceDeformation::interpolate(Implicit *imp)	
{
	if (imp==NULL)
		return;

	//[tk 5.19.05] I think everything is working except for this check to determing if *imp is of type RBF
	RBF* rbf_imp = dynamic_cast< RBF* >(imp);
	if (rbf_imp != NULL)
	//if (RBF* rbf_imp = dynamic_cast< RBF* >(imp))
	{
		//rbf_imp->interpolate(position->x, h, d_orient->n, flexible, position->changed);
		// [tk 5.19.05]
		// eventually there will need to be a vector of orientations to pass here to RBF to use as the normal constraints, but for now
		// they aren't used, so I'm passing the position vector again as the 3rd argument just to get the system working
		unsigned int hSz = h.size();
		
		// [tk 6.01.05]
		// I assume that particlePosition is not the only way to change the particle's position.  
		// This needs to be updated to be called after the particle system is updated in any way at all.  
		// Currently this code does not catch movements of the control points which should update through the GUI.
		// I need to track down what gets called when they are moved, and channel that signal to this point.
		if((position->changed) || (scalar->changed))
		{
			std::vector<double> interp_vals = scalar->getScalar();
			int sz = interp_vals.size();
			rbf_imp->interpolate(position->x, interp_vals, position->x, flexible, true);
			//rbf_imp->interpolate(position->x, scalar->getScalar(), position->x, flexible, position->changed);
			position->changed = false;
			scalar->changed = false;
		}
		return;
	}
	else
	{
	//continue;
	}

	int i,j,k;
	gmVector3 gradient;
	gmVector3 xi;
	int qlen = imp->qlen(); // # of parameters
	int n_free;              // # of free parameters
	int n = ps->size();      // # of control particles
	position->changed = false;


	// Check size of the flexible array
	if (flexible.size() != qlen)
	{
		flexible.resize(qlen);
		flexible = true;
		n_free = qlen;
	}
	else
	{
		n_free = 0;
		for (i = 0; i < qlen; i++)
		if (flexible[i]) 
			n_free++;
	}

	// Check for overconstrained case.
	if (n > n_free) 
		return;

	// Allocate some derivative vectors
	DoubleArray dFdQi(0.0, qlen);
	DoubleArray dFdQj(0.0, qlen);
	DoubleArray dqdt(0.0, qlen);

	/* Solve for Lagrangian multipliers
	* Equation (7) of WH.
	*
	* Sum_j (dFdQ[i].dFdQ[j]) lambda[j] = 
	*       dFdQ[i].dQdt + grad(x[i]).v[i] + phi*proc(x[i])
	*
	* Ax = b
	* x = vector of Lagrangians lambda[j], one per control particle
	* A = matrix of dot products of dFdQi with dFdQj
	* b = RHS of above equation
	*/

	TNT::Vector<double> b(n);
	TNT::Matrix<double> A(n,n);
	TNT::Vector<double> x(n);

	// Cycle through control particles
	for (i=0; i<n; i++)
	{
		// dFdQi = derivative of F wrt Q at control particle i
		imp->procq(position->getPosition(i), dFdQi);

		#ifdef DEBUG_MATRIX
		std::cerr << "dFdQ: ";
		for(int k=0;k<dFdQi.size();k++) 
		std::cerr << dFdQi[k] << " ";
		std::cerr << std::endl;
		#endif

		// Build element i of vector b
		b[i] = 0.0;

		// Add dFdQ[i].dQdt
		for (k = 0; k < qlen; k++)
			if (flexible[k])
				b[i] += dFdQi[k]*dqdt[k];

		// Add grad(x[i]).v[i] + phi*proc(x[i])
		b[i] += dot(imp->grad(position->getPosition(i)), velocity->v[i]) + phi * imp->proc(position->getPosition(i));

		// Build row i of matrix A
		for(j=0; j<n; j++)
		{
			A[i][j] = 0.0;

			// Get derivative of F wrt Q at particle j
			imp->procq(position->getPosition(j), dFdQj);

			// Only the subset of this array corresponding to flexible
			// parameters is used.
			for (k = 0; k < qlen; k++)
			if (flexible[k]) 
			A[i][j] += dFdQi[k]*dFdQj[k];
		} // end matrix loop
	} // end control particle loop

	#ifdef DEBUG_MATRIX
	std::cerr << A;
	#endif

	// Try solving Ax = b using Cholesky
	if (!solveCholesky(A, x, b))
	{
		// Let SVD take care of it
		if (!svd.solve(A, x, b)) 
			std::cerr << "cannot solving constraints!!" << std::endl;
	}

	#ifdef DEBUG_MATRIX
	TNT::Vector<double> r(n_controls);
	r = A*x - b;
	std::cerr << "residual: " << sqrt(TNT::dot_prod(r,r)) << std::endl;
	#endif

	// Update parameters
	for(j=0; j<n; j++)
	{
		imp->procq(position->getPosition(j), dFdQj);
		dFdQj *= (double)x[j];

		// Only modify those that are flexible
		for(k=0; k<qlen; k++) 
		if (flexible[k])
		dqdt[k] -= dFdQj[k];
	}

	// Apply the changes to the implicit surface parameter
	std::valarray<double> q(qlen);
	imp->getq(q);
	q += dqdt * (double)ps->dt;
	imp->setq(q);
}
#endif

/**
 * Solves a system of equations using Cholesky factorization.
 * @param A Matrix of coefficients.
 * @param x Vector of variables.
 * @param b Vector of values.
 * @return True if successful.
 */
bool SurfaceDeformation::solveCholesky(TNT::Matrix<double>& A,
									   TNT::Vector<double>& x,
                                       TNT::Vector<double>& b) 
{
	int n = x.size();
	TNT::Matrix<double> L(n, n);
	if (TNT::Cholesky_upper_factorization(A, L) != 0) return false;
	TNT::Vector<double> y = TNT::Lower_triangular_solve(L, b);
	TNT::Transpose_View<TNT::Matrix<double> > T = TNT::Transpose_view(L);
	x = TNT::Upper_triangular_solve(T, y);
	return true;
}
