/** @file CompactRBF.cpp 
 * Implementation of variational implicit surfaces using
 * compactly-supported radial basis functions.
 * Based off of the general RBF implementation by William Nagel.
 * @author Mike Flavin
 * @date April 13, 2004.
 */

//I've modified the KD-tree to support more general operations
//and have corrected several of the gradients, which were wrong
// -Scott Kircher

#include <algorithm>

#include "CompactRBF.h"

#include "tnt/tnt.h"
#include "tnt/cmat.h"
#include "tnt/vec.h"
#include "tnt/lu.h"
#include "Surface/gmTNTconvert.h"

REGISTER_IMPLICIT(CompactRBF,"Variational:CompactRBF");

/** Default empty constructor. This currently builds an RBF sphere consisting of six 0 centers
 * at (+/-1, +/-1, +/-1) and one -1 center at (0,0,0).
 */
CompactRBF::CompactRBF(void)
{
	constraints.resize(6);

	constraints[0] = RBFModelerConstraint(gmVector3(1.0,0.0,0.0), gmVector3(1.0,0.0,0.0));
	constraints[1] = RBFModelerConstraint(gmVector3(0.0,1.0,0.0), gmVector3(0.0,1.0,0.0));
	constraints[2] = RBFModelerConstraint(gmVector3(0.0,0.0,1.0), gmVector3(0.0,0.0,1.0));
	constraints[3] = RBFModelerConstraint(gmVector3(-1.0,0.0,0.0), gmVector3(-1.0,0.0,0.0));
	constraints[4] = RBFModelerConstraint(gmVector3(0.0,-1.0,0.0), gmVector3(0.0,-1.0,0.0));
	constraints[5] = RBFModelerConstraint(gmVector3(0.0,0.0,-1.0), gmVector3(0.0,0.0,-1.0));

	continuity = 2;
	support_radius = 1.5;
	param_changed = true;

	updateRBF();
}

/** Returns the number of parameters of the Compact RBF (which is 2).
 */
unsigned int CompactRBF::qlen()
{
	return 2;
}

/**
 * Loads parameters into q.
 * @param q An array of doubles, of size qlen, into which to copy the parameter array of the surface.
 */
void CompactRBF::getq(double *q)
{
	q[0] = static_cast<double>(continuity);
	q[1] = support_radius;
}

/** Stores parameters into q.
 * @param q An array of doubles, of size qlen, containing the values to store in the parameter vector.
 */
void CompactRBF::_setq(double *q)
{
	if(static_cast<int>(q[0])!=continuity) param_changed=true;
	if(support_radius!=q[1]) param_changed=true;
	continuity = static_cast<int>(q[0]);
	support_radius = q[1];
	
	if(param_changed)
	{
		updateRBF();
		param_changed=false;
	}
}

/**
 * Returns a list of parameter names.
 * @param qn an array of size qlen of strings listing names of parameters
 */
void CompactRBF::getqname(char **qn)
{
	qn[0] = "Continuity";
	qn[1] = "Support Radius";
}

/** 
 * Loads control particles into RBF vector and calls updateRBF().
 * @param positions A vector of constraint positions.
 * @param interp_vals A vector of desired function values for each constraint (determines interior, exterior, or surface constraint).
 * @param normals A vector of the desired normals for normal constraints.
 * @param flexible A vector which serves no discernable purpose for RBFs.
 * @param pos_changed Have any of the constraints changed since we last solved for a surface?  If not, do not re-solve.
 */
void CompactRBF::interpolate(std::vector<gmVector3> positions, std::vector<double> interp_vals, std::vector<gmVector3> normals, std::valarray<bool>& flexible, bool pos_changed)
{
	if(pos_changed || param_changed)
	{
		constraints.clear();
		constraints.resize(positions.size());

		for(unsigned int i = 0; i < positions.size(); i++)
		{
			// Normal constraint
			if(normals[i] != gmVector3(0.0, 0.0, 0.0))
				constraints[i] = RBFModelerConstraint(positions[i], normals[i], interp_vals[i], true);
			// Non-normal constraint
			else
				constraints[i] = RBFModelerConstraint(positions[i], normals[i], interp_vals[i], false);
		}
	
		updateRBF();
		param_changed = false;
	}
}

void CompactRBF::updateRBF(void)
{
  int i,j,n;

  convert_constraints();

  n = centers.size();

  // Insert all of the centers into the kd-tree.
  std::vector<gmVector3> points(n);
  for(int k = 0; k < n; k++)
  {
	  points[k]=centers[k].c;
  }
  support_tree.generate(points);

  TNT::Matrix<double> A(n + 4, n + 4);
  TNT::Vector<int> ipiv(n + 4);
  TNT::Vector<double> h(n + 4);

  for(i=0; i < n; i++) {
	for(j=0; j < n; j++) {
      // A[i][j] = (*this)[j].phi(convert((*this)[i].c));
      A[i][j] = phi(centers[j].c - centers[i].c);
	}

    A[n][i] = 1.0;
    A[i][n] = 1.0;

    for(j=0; j < 3; j++)
    {
      A[n + 1 + j][i] = centers[i].c[j];
      A[i][n + 1 + j] = centers[i].c[j];
    }
  }

  for(i = 0; i < n; i++)
    h[i] = centers[i].h;

  h[n] = 0.0;
  h[n + 1] = 0.0;
  h[n + 2] = 0.0;
  h[n + 3] = 0.0;

  TNT::LU_factor(A,ipiv);
  TNT::LU_solve(A,ipiv,h);

  for(i = 0; i < n; i++)
    centers[i].w = h[i];

  m_p[0] = h[n];
  m_p[1] = h[n+1];
  m_p[2] = h[n+2];
  m_p[3] = h[n+3];
}

/** Convert the "conceptual" constraints into actual RBF constraint points (with weights, etc).
 */      
void CompactRBF::convert_constraints()
{
	double normal_dist;

	// Adjust the distance between position and normal constraint point -
	//  they must always be within each other's support radius.
	if(support_radius >= 0.11)
		normal_dist = 0.1;
	else
		normal_dist = (0.9 * support_radius);

	centers.clear();

	// Convert all the constraint points.  Normal constraints become 2 separate "centers", with
	//  offset function value constraints.
	for(unsigned int i = 0; i < constraints.size(); i++)
	{
		centers.push_back(RBFPoint(constraints[i].c, constraints[i].h));
		if(constraints[i].nc)
		{
			centers.push_back(RBFPoint((constraints[i].c + (normal_dist * constraints[i].n)), (constraints[i].h + 1.0)));
		}
	}
}

#ifndef INTERVAL_EVAL_ONLY
/** Evaluate the variational implicit surface.
 * The surface consists of a list of radial basis functions and
 * a plane (the affine part).
 */
double CompactRBF::proc(const gmVector3 & x)
{
  double fVal = 0.0;
  
  std::vector<int> nz_indices;
  support_tree.radiusQuery(x, support_radius,nz_indices);
  
  // If the point is outside of the support radius of any point, make its value huge.
  //  This will cause ParticleValueDeath (if it is present) to destroy any floater
  //  outside of the support region, eliminating the "runaway floaters" problem.
  if(nz_indices.size() == 0)
	  return 100.0;

  for(unsigned int i = 0; i < nz_indices.size(); i++)
    fVal += centers[nz_indices[i]].w * phi(x - centers[nz_indices[i]].c);

  fVal += m_p[3] * x[2];
  fVal += m_p[2] * x[1];
  fVal += m_p[1] * x[0];
  fVal += m_p[0];

  return fVal;
}

/** Gradient of the sum of radial basis functions
 * phi = r^3
 * dphi = 3 r^2
 * Gradient is just magnitude r in the unitized direction of x.
 * grad phi = 3r x (since x = r x/||x||)
 */
gmVector3 CompactRBF::grad(const gmVector3 & x)
{
  gmVector3 g(0.0,0.0,0.0);

  std::vector<int> nz_indices;
  support_tree.radiusQuery(x, support_radius,nz_indices);

  for (unsigned int i = 0; i < nz_indices.size(); i++)
	g += centers[nz_indices[i]].w * gradphi(x - centers[nz_indices[i]].c);

  g += gmVector3(m_p[1],m_p[2],m_p[3]);

  return g;
}

gmMatrix3 CompactRBF::hess(const gmVector3 & x)
{
  gmMatrix3 g;

  std::vector<int> nz_indices;
  support_tree.radiusQuery(x, support_radius,nz_indices);

  for (unsigned int i = 0; i < nz_indices.size(); i++)
    g += centers[nz_indices[i]].w * hessphi(x - centers[nz_indices[i]].c);

  return g;
}
#endif

/** Evaluate the radial basis function.
 * For C^0 continuity, phi = (1 - ||x||)^2[+]
 * For C^2 continuity, phi = (1 - ||x||)^4[+] * (4 * ||x|| + 1)
 * For C^4 continuity, phi = (1 - ||x||)^6[+] * (35 * ||x||^2 + 18 * ||x|| + 3)
 * For C^6 continuity, phi = (1 - ||x||)^8[+] * (32 * ||x||^3 + 25 * ||x||^2 + 8 * ||x|| + 1)
 * @param x The vector from the location of the current constraint to the point at which we are evaluating the function.
 */
double CompactRBF::phi(const gmVector3 & x)
{
	if(x.lengthSquared() >= (support_radius * support_radius))
		return 0;

	if(continuity < 2)
		return phi_c0(x);
	else if(continuity < 4)
		return phi_c2(x);
	else if(continuity < 6)
		return phi_c4(x);
	else
		return phi_c6(x);
}

/** Evaluate the radial basis function over an interval.
 * @param b The interval over with to evaluate the function.
 */
/*Removed because I'm not sure they're correct, and we don't really need them --SIK
Intervald CompactRBF::phi(const Box<double>&  b)
{
	if(continuity < 2)
		return phi_c0(b);
	else if(continuity < 4)
		return phi_c2(b);
	else if(continuity < 6)
		return phi_c4(b);
	else
		return phi_c6(b);
}
*/

/** Evaluate the gradient of the RBF at a specific point in space.
 * @param x The vector from the location of the current constraint to the point at which we are evaluating the gradient.
 */
gmVector3 CompactRBF::gradphi(const gmVector3 & x)
{
if(x.lengthSquared() >= (support_radius * support_radius))
		return gmVector3(0.0, 0.0, 0.0);

	if(continuity < 2)
		return gradphi_c0(x);
	else if(continuity < 4)
		return gradphi_c2(x);
	else if(continuity < 6)
		return gradphi_c4(x);
	else
		return gradphi_c6(x);
}

/** Evaluate the gradient of the RBF over an interval.
 * @param b The interval over with to evaluate the gradient.
 */
/*Removed because I'm not sure they're correct, and we don't really need them --SIK
Box3d CompactRBF::gradphi(const Box<double>&  b)
{
	if(continuity < 2)
		return gradphi_c0(b);
	else if(continuity < 4)
		return gradphi_c2(b);
	else if(continuity < 6)
		return gradphi_c4(b);
	else
		return gradphi_c6(b);
}
*/

/** Evaluate the Hessian of the RBF (second-derivative matrix - roughly represents curvature of the surface at a point.)
 * @param x The vector from the location of the current constraint to the point at which we are evaluating the Hessian.
 */
gmMatrix3 CompactRBF::hessphi(const gmVector3 & x)
{
	if(x.lengthSquared() >= (support_radius * support_radius))
		return gmMatrix3();

	if(continuity < 2)
		return hessphi_c0(x);
	else if(continuity < 4)
		return hessphi_c2(x);
	else if(continuity < 6)
		return hessphi_c4(x);
	else
		return hessphi_c6(x);
}

/** Evaluate the Hessian of the radial basis function.
 *  The Hessian has a neat property - it can be broken into 2 parts:
 *   a uniform part that is the same for all entries (a fancy scalar
 *   times the outer product of the vector x with itself), and another
 *   part that is added to the partial derivative only along the diagonal.
 *   In addition, the second (non-uniform) part is just the identity matrix
 *   multiplied by the scalar length of the gradient.  This falls out
 *   of various derivative definitions, I'm sure.
 *
 * @param b The interval over with to evaluate the Hessian.
 */
/*Removed because I'm not sure they're correct, and we don't really need them --SIK
IMatrix3d CompactRBF::hessphi(const Box<double>&  b)
{ 
	if(continuity < 2)
		return hessphi_c0(b);
	else if(continuity < 4)
		return hessphi_c2(b);
	else if(continuity < 6)
		return hessphi_c4(b);
	else
		return hessphi_c6(b);
}
*/

/** The version of the RBF phi for C0 continuity.
 * @param x The vector from the location of the current constraint to the point at which we are evaluating the function.
 */
double CompactRBF::phi_c0(const gmVector3 & x)
{
	double rad = x.length() / support_radius;

	return (1.0 - rad) * (1.0 - rad);
}

/** The version of the RBF phi for C2 continuity.
 * @param x The vector from the location of the current constraint to the point at which we are evaluating the function.
 */
double CompactRBF::phi_c2(const gmVector3 & x)
{
	double rad = x.length() / support_radius;

	double temp= 1.0-rad;
	double term1 = temp*temp*temp*temp;
	double term2 = ((4.0 * rad) + 1.0);
	return (term1 * term2);
}

/** The version of the RBF phi for C4 continuity.
 * @param x The vector from the location of the current constraint to the point at which we are evaluating the function.
 */
double CompactRBF::phi_c4(const gmVector3 & x)
{
	double rad = x.length() / support_radius;

	double temp = 1.0-rad;
	double term1 = temp*temp*temp*temp*temp*temp;
	double term2 = ((35.0 * rad * rad) + (18.0 * rad) + 3.0);
	return (term1 * term2);
}

/** The version of the RBF phi for C6 continuity.
 * @param x The vector from the location of the current constraint to the point at which we are evaluating the function.
 */
double CompactRBF::phi_c6(const gmVector3 & x)
{
	double rad = x.length() / support_radius;

	double term1 = pow((1.0 - rad), 8.0);
	double term2 = ((32.0 * rad * rad * rad) + (25.0 * rad * rad) + (8.0 * rad) + 1.0);
	return (term1 * term2);
}

//The interval versions have been removed because I'm not sure their gradients/hessians are
//correct, and we don't really need the interval versions anyway --SIK
#if 0

/** The interval version of the RBF phi for C0 continuity.
 * @param b The interval over which to evaluate the function.
 */
Intervald CompactRBF::phi_c0(const Box<double>&  b)
{
	 Box3d b1 = Box3d(b[0],b[1],b[2]);
	 Intervald scaled = b1.length() / support_radius;
	 Intervald oneminusrad = (1.0 - scaled);
	 if(oneminusrad.low() < 0.0)
		 oneminusrad = Intervald(0.0, oneminusrad.high());
	 return oneminusrad.pow(2);
}

/** The interval version of the RBF phi for C2 continuity.
 * @param b The interval over which to evaluate the function.
 */
Intervald CompactRBF::phi_c2(const Box<double>&  b)
{
	Box3d b1 = Box3d(b[0],b[1],b[2]);
	Intervald scaled = b1.length() / support_radius;
	Intervald oneminusrad = (1.0 - scaled);
	if(oneminusrad.low() < 0.0)
		oneminusrad = Intervald(0.0, oneminusrad.high());
	Intervald temp1 = oneminusrad.pow(4);
	Intervald temp2 = (4.0 * scaled) + 1.0;
	return (temp1 * temp2);
}

/** The interval version of the RBF phi for C4 continuity.
 * @param b The interval over which to evaluate the function.
 */
Intervald CompactRBF::phi_c4(const Box<double>&  b)
{
	Box3d b1 = Box3d(b[0],b[1],b[2]);
	Intervald scaled = b1.length() / support_radius;
	Intervald oneminusrad = (1.0 - scaled);
	if(oneminusrad.low() < 0.0)
		oneminusrad = Intervald(0.0, oneminusrad.high());
	Intervald temp1 = oneminusrad.pow(8);
	Intervald temp2 = (32.0 * scaled.pow(3)) + (25.0 * scaled.pow(2)) + (8.0 * scaled) + 3.0;
	return (temp1 * temp2);
}

/** The interval version of the RBF phi for C6 continuity.
 * @param b The interval over which to evaluate the function.
 */
Intervald CompactRBF::phi_c6(const Box<double>&  b)
{
	Box3d b1 = Box3d(b[0],b[1],b[2]);
	Intervald scaled = b1.length() / support_radius;
	Intervald oneminusrad = (1.0 - scaled);
	if(oneminusrad.low() < 0.0)
		oneminusrad = Intervald(0.0, oneminusrad.high());
	Intervald temp1 = oneminusrad.pow(6);
	Intervald temp2 = (35.0 * scaled * scaled) + (18.0 * scaled) + 3.0;
	return (temp1 * temp2);
}
#endif

/** The version of the gradient for C0 continuity.
 * @param x The vector from the location of the current constraint to the point at which we are evaluating the gradient.
 */
gmVector3 CompactRBF::gradphi_c0(const gmVector3 & x)
{
	double rad = x.length() / support_radius;

	if(gmIsZero(rad))
		return x;
	else
		return (2.0 * (rad - 1.0) / support_radius) * x/x.length(); 
}

/** The version of the gradient for C2 continuity.
 * @param x The vector from the location of the current constraint to the point at which we are evaluating the gradient.
 */
gmVector3 CompactRBF::gradphi_c2(const gmVector3 & x)
{
/* 	double rad = x.length() / support_radius;
	double rm1_3 = pow((rad - 1.0), 3.0);
	return ((20.0 * rm1_3) * x);*/

	//Corrected version of gradient (the gradient of length() was ignored before) --SIK
	double rad = x.length() / support_radius;
	if(gmIsZero(rad)) return x;
	double temp = (1.0-rad);
	return 4.0*(temp*temp*temp*temp - temp*temp*temp*(4.0*rad+1))*x/(x.length()*support_radius);
}

/** The version of the gradient for C4 continuity.
 * @param x The vector from the location of the current constraint to the point at which we are evaluating the gradient.
 */
gmVector3 CompactRBF::gradphi_c4(const gmVector3 & x)
{
/*	double rad = x.length() / support_radius;
	double rm1_5 = pow((rad - 1.0), 5.0);
	return ((56.0 * rm1_5 * (5.0 * rad + 1.0)) * x);*/

	//Corrected version of gradient (the gradient of length() was ignored before) --SIK
	double rad = x.length() / support_radius;
	if(gmIsZero(rad)) return x;
	double temp = (1.0-rad);
	return -(280.0*rad*rad + 56.0*rad) * temp*temp*temp*temp*temp * x/(x.length()*support_radius);
}

/** The version of the gradient for C6 continuity.
 * @param x The vector from the location of the current constraint to the point at which we are evaluating the gradient.
 */
gmVector3 CompactRBF::gradphi_c6(const gmVector3 & x)
{
	double rad = x.length() / support_radius;
	if(gmIsZero(rad)) return x;
	double rm1_7 = pow((rad - 1.0), 7.0);
	double term2 = (16.0 * rad * rad) + (7.0 * rad) + 1.0;
	return ((22.0 * rm1_7 * term2) * x)/(x.length()*support_radius); //corrected by SIK

}


//The interval versions have been removed because I'm not sure their gradients/hessians are
//correct, and we don't really need the interval versions anyway --SIK
#if 0

/** The interval version of the gradient for C0 continuity.
 * @param b The interval over which to evaluate the gradient.
 */
Box3d CompactRBF::gradphi_c0(Box<double>b)
{
	Box3d b1 = Box3d(b[0],b[1],b[2]);
	Intervald rad = b1.length() / support_radius;
	if(rad.high() > 1.0)
		rad = Intervald(rad.low(), 1.0);
	Intervald radminus1 = rad - 1.0;
	Intervald term = 2.0 * radminus1 / rad;
	return (term * b1);
}

/** The interval version of the gradient for C2 continuity.
 * @param b The interval over which to evaluate the gradient.
 */
Box3d CompactRBF::gradphi_c2(Box<double>b)
{
	Box3d b1 = Box3d(b[0],b[1],b[2]);
	Intervald rad = b1.length() / support_radius;
	if(rad.high() > 1.0)
		rad = Intervald(rad.low(), 1.0);
	Intervald radminus1 = rad - 1.0;
	Intervald term = 20.0 * radminus1.pow(3);
	return (term * b1);
}

/** The interval version of the gradient for C4 continuity.
 * @param b The interval over which to evaluate the gradient.
 */
Box3d CompactRBF::gradphi_c4(Box<double>b)
{
	Box3d b1 = Box3d(b[0],b[1],b[2]);
	Intervald rad = b1.length() / support_radius;
	if(rad.high() > 1.0)
		rad = Intervald(rad.low(), 1.0);
	Intervald radminus1 = rad - 1.0;
	Intervald term1 = 56.0 * radminus1.pow(5);
	Intervald term2 = (5.0 * rad) + 1.0;
	return (term1 * term2 * b1);
}

/** The interval version of the gradient for C6 continuity.
 * @param b The interval over which to evaluate the gradient.
 */
Box3d CompactRBF::gradphi_c6(Box<double>b)
{
	Box3d b1 = Box3d(b[0],b[1],b[2]);
	Intervald rad = b1.length() / support_radius;
	if(rad.high() > 1.0)
		rad = Intervald(rad.low(), 1.0);
	Intervald radminus1 = rad - 1.0;
	Intervald term1 = (22.0 * radminus1.pow(7));
	Intervald term2 = (16.0 * rad.pow(2)) + (7.0 * rad) + 1.0;
	return (term1 * term2 * b1);
}
#endif

/** The version of the Hessian for C0 continuity.
 * @param x The vector from the location of the current constraint to the point at which we are evaluating the Hessian.
 */
gmMatrix3 CompactRBF::hessphi_c0(const gmVector3 & x)
{
	double rad = x.length() / support_radius;
	
	// Same cheat to avoid divide by zero as in regular RBFs.
	//  If someone wants to make a more elegant solution at some point,
	//  feel free.  For now, however, this works.
	if(gmIsZero(rad))
		rad = 0.000001;
		
	double weight1 = 1.0 / pow(rad, 3.0);
	double weight2 = 1.0 - (1.0 / rad);
	return 2.0 * ((weight1 * outer(x, x)) + (weight2 * gmMatrix3::identity()));
}

/** The version of the Hessian for C2 continuity.
 * @param x The vector from the location of the current constraint to the point at which we are evaluating the Hessian.
 */
gmMatrix3 CompactRBF::hessphi_c2(const gmVector3 & x)
{
	double rad = x.length() / support_radius;

	// Same cheat to avoid divide by zero as in regular RBFs.
	//  If someone wants to make a more elegant solution at some point,
	//  feel free.  For now, however, this works.
	if(gmIsZero(rad))
		rad = 0.000001;

	double radminus1 = rad - 1.0;
	double weight1 = 60.0 * radminus1 * radminus1 / rad;
	double weight2 = 20.0 * pow(radminus1, 3.0);
	return (weight1 * outer(x, x)) + (weight2 * gmMatrix3::identity());
}

/** The version of the Hessian for C4 continuity.
 * @param x The vector from the location of the current constraint to the point at which we are evaluating the Hessian.
 */
gmMatrix3 CompactRBF::hessphi_c4(const gmVector3 & x)
{
	double rad = x.length() / support_radius;
	double radminus1 = rad - 1.0;
	double weight1 = 1680.0 * pow(radminus1, 4.0);
	double weight2 = 56.0 * pow(radminus1, 5.0) * ((5.0 * rad) + 1.0);
	return (weight1 * outer(x, x)) + (weight2 * gmMatrix3::identity());
}

/** The version of the Hessian for C6 continuity.
 * @param x The vector from the location of the current constraint to the point at which we are evaluating the Hessian.
 */
gmMatrix3 CompactRBF::hessphi_c6(const gmVector3 & x)
{
	double rad = x.length() / support_radius;
	double radminus1 = rad - 1.0;
	double weight1 = 528.0 * pow(radminus1, 6.0) * (6.0 * rad + 1.0);
	double weight2 = 22.0 * pow(radminus1, 7.0) * ((16.0 * rad * rad) + (7.0 * rad) + 1.0);
	return (weight1 * outer(x, x)) + (weight2 * gmMatrix3::identity());
}

//The interval versions have been removed because I'm not sure their gradients/hessians are
//correct, and we don't really need the interval versions anyway --SIK
#if 0

/** The interval version of the Hessian for C0 continuity.
 * @param b The interval over which to evaluate the Hessian.
 */
IMatrix3d CompactRBF::hessphi_c0(const Box<double>&  b)
{
	Box3d b1 = Box3d(b[0],b[1],b[2]);
	Intervald rad = b1.length() / support_radius;
	if(rad.high() > 1.0)
		rad = Intervald(rad.low(), 1.0);
	IMatrix3d uniform;
	IMatrix3d diag;
	uniform = (2.0 / rad.pow(3.0)) * IMatrix(outer(b1, b1));
	diag[0][0] = diag[1][1] = diag[2][2] = 2.0 - (2.0 / rad);
	return uniform + diag;
}

/** The interval version of the Hessian for C2 continuity.
 * @param b The interval over which to evaluate the Hessian.
 */
IMatrix3d CompactRBF::hessphi_c2(const Box<double>&  b)
{
	Box3d b1 = Box3d(b[0],b[1],b[2]);
	Intervald rad = b1.length() / support_radius;
	if(rad.high() > 1.0)
		rad = Intervald(rad.low(), 1.0);
	Intervald radminus1 = rad - 1.0;
	IMatrix3d uniform;
	IMatrix3d diag;
	uniform = (60.0 * radminus1.pow(2) / rad) * IMatrix(outer(b1, b1));
	diag[0][0] = diag[1][1] = diag[2][2] = 20.0 * radminus1.pow(3);
	return uniform + diag;
}

/** The interval version of the Hessian for C4 continuity.
 * @param b The interval over which to evaluate the Hessian.
 */
IMatrix3d CompactRBF::hessphi_c4(const Box<double>&  b)
{
	Box3d b1 = Box3d(b[0],b[1],b[2]);
	Intervald rad = b1.length() / support_radius;
	if(rad.high() > 1.0)
		rad = Intervald(rad.low(), 1.0);
	Intervald radminus1 = rad - 1.0;
	IMatrix3d uniform;
	IMatrix3d diag;
	uniform = (1680.0 * radminus1.pow(4)) * IMatrix(outer(b1, b1));
	diag[0][0] = diag[1][1] = diag[2][2] = 56.0 * radminus1.pow(5) * (5.0 * rad + 1.0);
	return uniform + diag;
}

/** The interval version of the Hessian for C4 continuity.
 * @param b The interval over which to evaluate the Hessian.
 */
IMatrix3d CompactRBF::hessphi_c6(const Box<double>&  b)
{
	Box3d b1 = Box3d(b[0],b[1],b[2]);
	Intervald rad = b1.length() / support_radius;
	if(rad.high() > 1.0)
		rad = Intervald(rad.low(), 1.0);
	Intervald radminus1 = rad - 1.0;
	IMatrix3d uniform;
	IMatrix3d diag;
	uniform = (528.0 * radminus1.pow(6) * (6.0 * rad + 1.0)) * IMatrix(outer(b1, b1));
	diag[0][0] = diag[1][1] = diag[2][2] = 22.0 * radminus1.pow(7) * ((16.0 * rad.pow(2)) + (7.0 * rad) + 1.0);
	return uniform + diag;
}
#endif
