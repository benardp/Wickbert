
/**
* Implementation of KDPropagation.
* @file KDPropagation.cpp
* @author John Hart
*/

#include "KDPropagation.h"
// #include "Implicit/Variational/RBF.h"
#include "ParticlePosition.h"
#include "ParticleVelocity.h"

#include <fstream>
#include <sstream>
using namespace std;

// Witkin-Heckbert recommended default values
#define WH_PHI 15.0

REGISTER_PARTICLESTUFF(KDPropagation,"Behavior:KDPropagation");

//void KDPropagation::attachAttributes()
//{
//	attachAttribute(imp_int,std::string("ImplicitInterrogator"));
//	attachAttribute(p_orient,std::string("ParticleOrientation"));
//}

/**
* Creates a surface Propagation behavior for Particles p.
* @param ps The owning particle system.
*
* Sets phi (the feedback constant) to 15.
* 
*/
KDPropagation::KDPropagation(Particles *ps, std::string name)
: SurfacePropagation(ps, std::string("KDPropagation"))
{
	/*new PSParamDouble(this,&speedfac,0.0,"speedfac","Speed Factor",
		"Global speed factor of all speed functions. "
		"Set to zero to disable growth.");
	new PSParamDouble(this,&constfac,1.0,"constfac","Constant Speed Term",
		"Term controlling constant growth rate.");
	new PSParamDouble(this,&gaussfac,0.0,"gaussfac","Gaussian Curvature Speed Term",
		"Term controlling proportion of growth rate due to Gassian curvature.");
	new PSParamDouble(this,&meanfac,0.0,"meanfac","Mean Curvature Speed Term",
		"Term controlling proportion of growth rate due to mean curvature.");

	new Attached<ImplicitInterrogator>(this,&imp_int);
	new Attached<ParticleOrientation>(this,&p_orient);
	new Attached<ParticlePosition>(this,&position);
	new Attached<ParticleVelocity>(this,&velocity);*/

	new PSParamString(this, &ptcloud_file, "", "", "", "input", "input", "input");
	new PSParamBool(this, &reload, false, "reload", "reload", "reload");

	loaded = false;
	//load_point_cloud();
}

void KDPropagation::load_point_cloud()
{
	//	load the point cloud
	//ifstream fs("c:\\xinlai\\inactives\\meshlib\\v2.pcl");
	ifstream fs(ptcloud_file.c_str());
	istringstream str;
	if( !fs )
		loaded = false;
	string line;
	double f;
	kdpointt pt;
	while( getline(fs, line) ) {
		str.clear();
		str.str(line);
		for(int i = 0; i < 3; i++) {
			str >> pt.pos[i];
		}
		str >> f;
		for(int i = 0; i < 3; i++) {
			str >> pt.norm[i];
		}
		points.push_back(pt);
	}
	//	construct the kd tree
	kd.construct(points);
	reload = false;
	loaded = true;
}

double KDPropagation::speed(int i)
{
	double s;

	s = constfac;

	if (!imp_int)
		return s;
	Implicit *imp = imp_int->getImplicit();
	if (!imp)
		return s;
	
//	s += gaussfac*imp->gaussianCurvature(position->x[i]);
//	s += meanfac*imp->meanCurvature(position->x[i]);

	//s *= speedfac;
	if( reload ) 
		load_point_cloud();
	if( loaded ) {
		double query[3];
		query[0] = position->x[i][0];
		query[1] = position->x[i][1];
		query[2] = position->x[i][2];
		int id = kd.search_nearest(query, points);
		double v[3];
		v[0] = points[id].pos[0] - query[0];
		v[1] = points[id].pos[1] - query[1];
		v[2] = points[id].pos[2] - query[2];
		s = v[0] * points[id].norm[0] + v[1] * points[id].norm[1] + v[2] * points[id].norm[2];
	}
	return s;
}


//void KDPropagation::applyConstraintOLD()
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


//// In the non-trivial case where there are more particles than parameters, we have
//// an overdetermined system.  We solve it with the SVD.  The old (unused) implementation
//// is in applyConstraintOLD
//void KDPropagation::applyConstraint()
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
//void KDPropagation::integrate()
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
//bool KDPropagation::solveCholesky(TNT::Matrix<double>& A,
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
//bool KDPropagation::solveLU(TNT::Matrix<double>& A, TNT::Vector<double>& b)
//{
//	TNT::Vector<TNT::Subscript> ipiv;
//	if ( TNT::LU_factor(A, ipiv) != 0 ) return false;
//	TNT::LU_solve(A, ipiv, b);
//
//	return true;
//}
