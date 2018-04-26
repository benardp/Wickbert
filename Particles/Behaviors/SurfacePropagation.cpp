/**
* Implementation of SurfacePropagation.
* @file SurfacePropagation.cpp
* @author John Hart
*/

#include "SurfacePropagation.h"
// #include "Implicit/Variational/RBF.h"
#include "ParticlePosition.h"
#include "ParticleVelocity.h"
#include "ParticleScalar.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>

// Witkin-Heckbert recommended default values
#define WH_PHI 15.0

REGISTER_PARTICLESTUFF(SurfacePropagation,"Behavior:SurfacePropagation");

/**
* Creates a surface Propagation behavior for Particles p.
* @param ps The owning particle system.
*
* Sets phi (the feedback constant) to 15.
* 
*/
SurfacePropagation::SurfacePropagation(Particles *ps, std::string name)
: ParticleBehavior(ps, std::string("SurfacePropagation"))
{
	new PSParamDouble(this,&speedfac,0.0,"speedfac","Speed Factor",
		"Global speed factor of all speed functions. "
		"Set to zero to disable growth.");
	new PSParamDouble(this,&constfac,1.0,"constfac","Constant Speed Term",
		"Term controlling constant growth rate.");
#if 0
	gaussfac = meanfac = 0.0;
#else
	new PSParamDouble(this,&gaussfac,0.0,"gaussfac","Gaussian Curvature Speed Term",
		"Term controlling proportion of growth rate due to Gassian curvature.");
	new PSParamDouble(this,&meanfac,0.0,"meanfac","Mean Curvature Speed Term",
		"Term controlling proportion of growth rate due to mean curvature.");
#endif
	new PSParamDouble(this,&tolerance,0.00001,"tolerance","Integration Tolerance","Term controlling the accuracy of the numerical integration");

	new Attached<ImplicitInterrogator>(this,&imp_int);	
	new Attached<ImplicitInterrogator>(this,&speed_int,"SpeedInterrogator","speed","Speed Func.",
		"Interrogator for an implicit used to control propagation speed.");
	new Attached<ParticleNormal>(this,&p_orient);
	new Attached<ParticlePosition>(this,&position);
	new Attached<ParticleVelocity>(this,&velocity);
	new Attached<ParticleScalar>(this,&perr);

	qdotPrev = 0;
}

double SurfacePropagation::speed(int i)
{
	double s;

	s = constfac;

#if 0
	if (!imp_int)
		return s;
	Implicit *imp = imp_int->getImplicit();
	if (!imp)
		return s;
	
	s += gaussfac*imp->gaussianCurvature(position->x[i]);
	s += meanfac*imp->meanCurvature(position->x[i]);
#else
	Implicit *simp = speed_int->getImplicit();
	if (simp && position)
		s += simp->proc(position->x[i]);
#endif

//	s *= speedfac;
	return s;
}


void SurfacePropagation::applyConstraintOLD()
{
	if (!imp_int) return;

	// imp_int->getImplicit()->interpolate(ps,flexible);
	Implicit *imp = imp_int->getImplicit();
	if (!imp) return;

	if (speedfac == 0.0)
		return;
	
	int i,j,k;
	gmVector3 gradient;
	gmVector3 xi;
	int qlen = imp->qlen(); // # of parameters
	int n = ps->size();      // # of control particles

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
	for (i=0; i<n; i++) {
		// dFdQi = derivative of F wrt Q at control particle i
		imp->procq(position->x[i], dFdQi);
	    
#ifdef DEBUG_MATRIX
		std::cerr << "dFdQ: ";
		for(int k=0;k<dFdQi.size();k++) 
			std::cerr << dFdQi[k] << " ";
		std::cerr << std::endl;
#endif
      
		// b = - speedfac * speed function * grad(x[i]) magnitude * speed function
		b[i] = -speedfac*speed(i)*imp->grad(position->x[i]).length();

		// Build row i of matrix A
		for(j=0; j<n; j++) {
			A[i][j] = 0.0;
	          
			// Get derivative of F wrt Q at particle j
			imp->procq(position->x[j], dFdQj);

			// Aij = dFdQi . dFdQj
			for (k = 0; k < qlen; k++)
				A[i][j] += dFdQi[k]*dFdQj[k];
	    } // end matrix loop
    } // end control particle loop
  
#ifdef DEBUG_MATRIX
	std::cerr << A;
#endif
  
	// Try solving Ax = b using Cholesky
	if (!solveCholesky(A, x, b)) {
		//std::cerr << "Using SVD... " << std::endl;
		// Let SVD take care of it
		SVD svd;
		if (!svd.solve(A, x, b)) 
			std::cerr << "problem solving constraints!!" << std::endl;
	}
  
#ifdef DEBUG_MATRIX
	TNT::Vector<double> r(n);
	r = A*x - b;
	double residual = sqrt(TNT::dot_prod(r, r));
	std::cerr << "residual: " << sqrt(TNT::dot_prod(r,r)) << std::endl;
#endif
  
	// Update parameters
	for(j=0; j<n; j++) {
		imp->procq(position->x[j], dFdQj);
		dFdQj *= (double)x[j];
	      
		for(k=0; k<qlen; k++) 
			dqdt[k] += dFdQj[k];
    }

	// Apply the changes to the implicit surface parameter
	std::valarray<double> q(qlen);
	imp->getq(q);
	q += dqdt * (double)ps->dt;
	imp->setq(q);
}


// In the non-trivial case where there are more particles than parameters, we have
// an overdetermined system.  We solve it with the SVD.  The old (unused) implementation
// is in applyConstraintOLD
#if 0
void SurfacePropagation::applyConstraint()
{
	if (!imp_int) return;

	Implicit *imp = imp_int->getImplicit();
	if (!imp) return;

	if (speedfac == 0.0)
		return;
	
	int i,j,k;
	int qlen = imp->qlen(); // # of parameters
	int n = ps->size();      // # of control particles

	// Allocate some derivative vectors
	DoubleArray dqdt(0.0, qlen);

	TNT::Vector<double> b(n);
	TNT::Matrix<double> A(n, qlen);
	TNT::Vector<double> x(qlen);

	// Cycle through control particles
	for (i=0; i<n; i++) {
		// dFdQi = derivative of F wrt Q at control particle i
//		imp->procq(position->x[i], dFdQi);
		imp->procq(position->x[i], A[i]);

		b[i] = -speedfac*speed(i)*imp->grad(position->x[i]).length();
	}

	
	// Use the SVD to solve the least squares problem
	SVD svd;
	if ( !svd.solve(A, x, b) )
		std::cout << "SVD solve failed." << std::endl;
	double a[10];
	for (j=0; j < qlen; j++) {
		a[j] = (double)x[j];
		dqdt[j] = (double)x[j];
	}

	// Apply the changes to the implicit surface parameter
	std::valarray<double> q(qlen);
	imp->getq(q);
	q += dqdt * (double)ps->dt;
	imp->setq(q);
}
#endif

#if 1
void SurfacePropagation::applyConstraint()
{
	if (!imp_int) return;

	Implicit *imp = imp_int->getImplicit();
	if (!imp) return;

	if (speedfac == 0.0)
		return;
	
	unsigned int i,j;
	unsigned int qlen = imp->qlen(); // # of parameters
	unsigned int n = ps->size();      // # of control particles

	// qdotPrev is the previous qdot.  We use it for a backwards difference
	// approximation to q'' in order to ensure stability and accuracy
	bool qPrevValid(true);
	if ( !qdotPrev || qdotPrev->size() != qlen )
	{
		// reset qdotPrev
		delete qdotPrev;
		qdotPrev = new DoubleArray(0.0, qlen);
		qPrevValid = false;
	}


	// Allocate some derivative vectors
	DoubleArray dqdt(0.0, qlen);

	double* A = new double[n * qlen];
	double* b = new double[n];
	for(i = 0; i < n; i++) {
		imp->procq(position->x[i], A + i * qlen);
		b[i] = -speedfac*speed(i)*imp->grad(position->x[i]).length();
		gmVector3 m = speedfac*speed(i)*imp->grad(position->x[i]);
		m.normalize();
		gmVector3 mm = velocity->v[i];
		mm.normalize();
		if ( dot(mm, m) < 0.90 )
			perr->setScalar(i, 0.5);
		else
			perr->setScalar(i, 1.0);
//		perr->setScalar(i, dot(mm, m));
	}
#if 0
	std::ofstream out("c:\\xinlai\\pointclouds\\A.txt");
	std::ofstream out2("c:\\xinlai\\pointclouds\\b.txt");
	for(unsignedint i = 0; i < n; i++)
	{
		for(int j = 0; j < qlen; j++)
			out << A[i * qlen + j] << " ";
		out << std::endl;
		out2 << b[i] << std::endl;
	}
#endif
	gsl_matrix_view gsl_A = gsl_matrix_view_array(A, n, qlen);
	gsl_matrix* V = gsl_matrix_alloc(qlen, qlen);
	gsl_matrix* work_X = gsl_matrix_alloc(qlen, qlen);
	gsl_vector* S = gsl_vector_alloc(qlen);
	gsl_vector* work = gsl_vector_alloc(qlen);
	gsl_linalg_SV_decomp(&gsl_A.matrix, V, S, work);
	//gsl_linalg_SV_decomp_mod(&gsl_A.matrix, work_X, V, S, work);
	//gsl_linalg_SV_decomp_jacobi(&gsl_A.matrix, V, S);


	gsl_vector_view gsl_b = gsl_vector_view_array(b, n);
	gsl_vector* x = gsl_vector_alloc(qlen);
	for(unsigned int i = 0; i < qlen; i++) {
		x->data[i] = 0.0;
		if( S->data[i] < 0.01 )
			S->data[i] = 0.0;
	}

	gsl_linalg_SV_solve(&gsl_A.matrix, V, S, &gsl_b.vector, x);
	

#if 0
	double norm = 0.0;
	for (j=0; j < qlen; j++) {
		double dq = (double)gsl_vector_get(x, j);
	norm += dq * dq;
	}
	norm = sqrt(norm);
#endif
	for(j = 0; j < qlen; j++) {
		dqdt[j] = (double)gsl_vector_get(x, j);
	}

	gsl_vector_free(x);
	gsl_vector_free(S);
	gsl_vector_free(work);
	gsl_matrix_free(work_X);
	gsl_matrix_free(V);
	delete[] A;
	delete[] b;

	// Apply the changes to the implicit surface parameter
	std::valarray<double> q(qlen);
	imp->getq(q);
	
	double stepsize;
	if ( qPrevValid )
	{
		// backward difference approximation to q''
		double qddNorm(0.0);
		for (unsigned int i=0; i < qlen; i++) {
			double x = dqdt[i] - (*qdotPrev)[i];
			qddNorm += x * x;
		}
		qddNorm = sqrt(qddNorm) / ps->dt;
		stepsize = sqrt(2 * tolerance / qddNorm);
	} else
		stepsize = ps->dt;

	q += dqdt * stepsize;
	imp->setq(q);

	for (unsigned int i=0; i < qlen; i++)
		(*qdotPrev)[i] = dqdt[i];
}

#endif

void SurfacePropagation::integrate()
{
}

/**
* Solves a system of equations using Cholesky factorization.
* @param A Matrix of coefficients.
* @param b Vector of values.	x->data[16]	0.40654805973887242	double

* @param x Vector of variables.
* @returns True if successful.
*/
bool SurfacePropagation::solveCholesky(TNT::Matrix<double>& A,
									   TNT::Vector<double>& x,
									   TNT::Vector<double>& b) {
	
	int n = ps->size();
	
	TNT::Matrix<double> L(n, n);
	
	if (TNT::Cholesky_upper_factorization(A, L) != 0) return false;
	
	TNT::Vector<double> y = TNT::Lower_triangular_solve(L, b);
	
	TNT::Transpose_View<TNT::Matrix<double> > T = TNT::Transpose_view(L);
	
	x = TNT::Upper_triangular_solve(T, y);
	
	return true;
}


bool SurfacePropagation::solveLU(TNT::Matrix<double>& A, TNT::Vector<double>& b)
{
	TNT::Vector<TNT::Subscript> ipiv;
	if ( TNT::LU_factor(A, ipiv) != 0 ) return false;
	TNT::LU_solve(A, ipiv, b);

	return true;
}

//
///**
//* Implementation of SurfacePropagation.
//* @file SurfacePropagation.cpp
//* @author John Hart
//*/
//
//#include "SurfacePropagation.h"
//// #include "Implicit/Variational/RBF.h"
//#include "ParticlePosition.h"
//#include "ParticleVelocity.h"
//
//#include <fstream>
//#include <sstream>
//using namespace std;
//
//// Witkin-Heckbert recommended default values
//#define WH_PHI 15.0
//
//REGISTER_PARTICLESTUFF(SurfacePropagation,"Behavior:SurfacePropagation");
//
//
///**
//* Creates a surface Propagation behavior for Particles p.
//* @param ps The owning particle system.
//*
//* Sets phi (the feedback constant) to 15.
//* 
//*/
//SurfacePropagation::SurfacePropagation(Particles *ps, std::string name)
//: ParticleBehavior(ps, std::string("SurfacePropagation"))
//{
//	new PSParamDouble(this,&speedfac,0.0,"speedfac","Speed Factor",
//		"Global speed factor of all speed functions. "
//		"Set to zero to disable growth.");
//	new PSParamDouble(this,&constfac,1.0,"constfac","Constant Speed Term",
//		"Term controlling constant growth rate.");
//	new PSParamDouble(this,&gaussfac,0.0,"gaussfac","Gaussian Curvature Speed Term",
//		"Term controlling proportion of growth rate due to Gassian curvature.");
//	new PSParamDouble(this,&meanfac,0.0,"meanfac","Mean Curvature Speed Term",
//		"Term controlling proportion of growth rate due to mean curvature.");
//
//	new Attached<ImplicitInterrogator>(this,&imp_int);
//	new Attached<ParticleOrientation>(this,&p_orient);
//	new Attached<ParticlePosition>(this,&position);
//	new Attached<ParticleVelocity>(this,&velocity);
//
//	//	load the point cloud
//	ifstream fs("c:\\xinlai\\inactives\\meshlib\\v2.pcl");
//	istringstream str;
//	if( !fs )
//		exit(1);
//	string line;
//	double f;
//	kdpointt pt;
//	while( getline(fs, line) ) {
//		str.clear();
//		str.str(line);
//		for(int i = 0; i < 3; i++) {
//			str >> pt.pos[i];
//		}
//		str >> f;
//		for(int i = 0; i < 3; i++) {
//			str >> pt.norm[i];
//		}
//		points.push_back(pt);
//	}
//	//	construct the kd tree
//	kd.construct(points);
//}
//
//double SurfacePropagation::speed(int i)
//{
//	double s;
//
//	s = constfac;
//
//	if (!imp_int)
//		return s;
//	Implicit *imp = imp_int->getImplicit();
//	if (!imp)
//		return s;
//	
////	s += gaussfac*imp->gaussianCurvature(position->x[i]);
////	s += meanfac*imp->meanCurvature(position->x[i]);
//
//	//s *= speedfac;
//	double query[3];
//	query[0] = position->x[i][0];
//	query[1] = position->x[i][1];
//	query[2] = position->x[i][2];
//	int id = kd.search_nearest(query, points);
//	double v[3];
//	v[0] = points[id].pos[0] - query[0];
//	v[1] = points[id].pos[1] - query[1];
//	v[2] = points[id].pos[2] - query[2];
//	s = v[0] * points[id].norm[0] + v[1] * points[id].norm[1] + v[2] * points[id].norm[2];
//	return s;
//}
//
//
//void SurfacePropagation::applyConstraintOLD()
//{
//	if (!imp_int) return;
//
//	// imp_int->getImplicit()->interpolate(ps,flexible);
//	Implicit *imp = imp_int->getImplicit();
//	if (!imp) return;
//
//	if (speedfac == 0.0)
//		return;
//	
//	int i,j,k;
//	gmVector3 gradient;
//	gmVector3 xi;
//	int qlen = imp->qlen(); // # of parameters
//	int n = ps->size();      // # of control particles
//	double phi = 15.0;
//
//
//	// Allocate some derivative vectors
//	DoubleArray dFdQi(0.0, qlen);
//	DoubleArray dFdQj(0.0, qlen);
//	DoubleArray dqdt(0.0, qlen);
//  
//	/* Solve for Lagrangian multipliers
//	* Equation (7) of WH.
//	*
//	* Sum_j (dFdQ[i].dFdQ[j]) lambda[j] = 
//	*       dFdQ[i].dQdt + grad(x[i]).v[i] + phi*proc(x[i])
//	*
//	* Ax = b
//	* x = vector of Lagrangians lambda[j], one per control particle
//	* A = matrix of dot products of dFdQi with dFdQj
//	* b = RHS of above equation
//	*/
//    
//	TNT::Vector<double> b(n);
//	TNT::Matrix<double> A(n,n);
//	TNT::Vector<double> x(n);
//
//	
//	// Cycle through control particles
//	for (i=0; i<n; i++) {
//		// dFdQi = derivative of F wrt Q at control particle i
//		imp->procq(position->x[i], dFdQi);
//	    
//#ifdef DEBUG_MATRIX
//		std::cerr << "dFdQ: ";
//		for(int k=0;k<dFdQi.size();k++) 
//			std::cerr << dFdQi[k] << " ";
//		std::cerr << std::endl;
//#endif
//      
//		// b = - speedfac * speed function * grad(x[i]) magnitude * speed function
//		b[i] = -speedfac*speed(i)*imp->grad(position->x[i]).length();
//
//		// Build row i of matrix A
//		for(j=0; j<n; j++) {
//			A[i][j] = 0.0;
//	          
//			// Get derivative of F wrt Q at particle j
//			imp->procq(position->x[j], dFdQj);
//
//			// Aij = dFdQi . dFdQj
//			for (k = 0; k < qlen; k++)
//				A[i][j] += dFdQi[k]*dFdQj[k];
//	    } // end matrix loop
//    } // end control particle loop
//  
//#ifdef DEBUG_MATRIX
//	std::cerr << A;
//#endif
//  
//	// Try solving Ax = b using Cholesky
//	if (!solveCholesky(A, x, b)) {
//		//std::cerr << "Using SVD... " << std::endl;
//		// Let SVD take care of it
//		SVD svd;
//		if (!svd.solve(A, x, b)) 
//			std::cerr << "problem solving constraints!!" << std::endl;
//	}
//  
//#ifdef DEBUG_MATRIX
//	TNT::Vector<double> r(n_controls);
//	r = A*x - b;
//	std::cerr << "residual: " << sqrt(TNT::dot_prod(r,r)) << std::endl;
//#endif
//  
//	// Update parameters
//	for(j=0; j<n; j++) {
//		imp->procq(position->x[j], dFdQj);
//		dFdQj *= (double)x[j];
//	      
//		for(k=0; k<qlen; k++) 
//			dqdt[k] += dFdQj[k];
//    }
//
//	// Apply the changes to the implicit surface parameter
//	std::valarray<double> q(qlen);
//	imp->getq(q);
//	q += dqdt * (double)ps->dt;
//	imp->setq(q);
//}
//
//
//// In the non-trivial case where there are more particles than parameters, we have
//// an overdetermined system.  We solve it with the SVD.  The old (unused) implementation
//// is in applyConstraintOLD
//void SurfacePropagation::applyConstraint()
//{
//	if (!imp_int) return;
//
//	Implicit *imp = imp_int->getImplicit();
//	if (!imp) return;
//
//	if (speedfac == 0.0)
//		return;
//	
//	int i,j,k;
//	int qlen = imp->qlen(); // # of parameters
//	int n = ps->size();      // # of control particles
//
//	// Allocate some derivative vectors
//	DoubleArray dqdt(0.0, qlen);
//
//	TNT::Vector<double> b(n);
//	TNT::Matrix<double> A(n, qlen);
//	TNT::Vector<double> x(qlen);
//
//	// Cycle through control particles
//	for (i=0; i<n; i++) {
//		// dFdQi = derivative of F wrt Q at control particle i
////		imp->procq(position->x[i], dFdQi);
//		imp->procq(position->x[i], A[i]);
//
//		b[i] = -speedfac*speed(i)*imp->grad(position->x[i]).length();
//	}
//
//	
//	// Use the SVD to solve the least squares problem
//	SVD svd;
//	if ( !svd.solve(A, x, b) )
//		std::cout << "SVD solve failed." << std::endl;
//
//	for (j=0; j < qlen; j++)
//		dqdt[j] = (double)x[j];
//
//	// Apply the changes to the implicit surface parameter
//	std::valarray<double> q(qlen);
//	imp->getq(q);
//	q += dqdt * (double)ps->dt;
//	imp->setq(q);
//}
//
//
//
//void SurfacePropagation::integrate()
//{
//}
//
///**
//* Solves a system of equations using Cholesky factorization.
//* @param A Matrix of coefficients.
//* @param b Vector of values.
//* @param x Vector of variables.
//* @returns True if successful.
//*/
//bool SurfacePropagation::solveCholesky(TNT::Matrix<double>& A,
//									   TNT::Vector<double>& x,
//									   TNT::Vector<double>& b) {
//	
//	int n = ps->size();
//	
//	TNT::Matrix<double> L(n, n);
//	
//	if (TNT::Cholesky_upper_factorization(A, L) != 0) return false;
//	
//	TNT::Vector<double> y = TNT::Lower_triangular_solve(L, b);
//	
//	TNT::Transpose_View<TNT::Matrix<double> > T = TNT::Transpose_view(L);
//	
//	x = TNT::Upper_triangular_solve(T, y);
//	
//	return true;
//}
//
//
//bool SurfacePropagation::solveLU(TNT::Matrix<double>& A, TNT::Vector<double>& b)
//{
//	TNT::Vector<TNT::Subscript> ipiv;
//	if ( TNT::LU_factor(A, ipiv) != 0 ) return false;
//	TNT::LU_solve(A, ipiv, b);
//
//	return true;
//}
