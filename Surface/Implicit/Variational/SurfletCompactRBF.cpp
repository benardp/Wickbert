/** @file SurfletCompactRBF.cpp 
 * Implementation of variational implicit surfaces using
 * surflet-style compactly-supported radial basis functions.
 * Extension of Mike Flavin's implementation of Compact RBFs.
 * @author Scott Kircher
 * @date November 7, 2004.
 */

#include <algorithm>

#include "SurfletCompactRBF.h"

#include "tnt/tnt.h"
#include "tnt/cmat.h"
#include "tnt/vec.h"
#include "tnt/lu.h"
#include "tnt/jama_qr.h"
#include "Surface/gmTNTconvert.h"

/*************************SURFLETCOMPACTRBF IMPLEMENTATION*************/

REGISTER_IMPLICIT(SurfletCompactRBF,"Variational:SurfletCompactRBF");

/** Default empty constructor. This currently builds an RBF sphere consisting of twelve 0 centers
 * at (+/-1, +/-1, +/-1)
 */
SurfletCompactRBF::SurfletCompactRBF()
{
	constraints.resize(12);

	constraints[0] = RBFModelerConstraint(gmVector3(1.0,0.0,0.0), gmVector3(1.0,0.0,0.0));
	constraints[1] = RBFModelerConstraint(gmVector3(0.0,1.0,0.0), gmVector3(0.0,1.0,0.0));
	constraints[2] = RBFModelerConstraint(gmVector3(0.0,0.0,1.0), gmVector3(0.0,0.0,1.0));
	constraints[3] = RBFModelerConstraint(gmVector3(-1.0,0.0,0.0), gmVector3(-1.0,0.0,0.0));
	constraints[4] = RBFModelerConstraint(gmVector3(0.0,-1.0,0.0), gmVector3(0.0,-1.0,0.0));
	constraints[5] = RBFModelerConstraint(gmVector3(0.0,0.0,-1.0), gmVector3(0.0,0.0,-1.0));
	constraints[6] = RBFModelerConstraint(gmVector3(1.0,1.0,0.0), gmVector3(1.0,1.0,0.0).normalize());
	constraints[7] = RBFModelerConstraint(gmVector3(0.0,1.0,1.0), gmVector3(0.0,1.0,1.0).normalize());
	constraints[8] = RBFModelerConstraint(gmVector3(-1.0,0.0,1.0), gmVector3(0.0,0.0,1.0).normalize());
	constraints[9] = RBFModelerConstraint(gmVector3(-1.0,-1.0,0.0), gmVector3(-1.0,-1.0,0.0).normalize());
	constraints[10] = RBFModelerConstraint(gmVector3(0.0,-1.0,-1.0), gmVector3(0.0,-1.0,-1.0).normalize());
	constraints[11] = RBFModelerConstraint(gmVector3(1.0,0.0,-1.0), gmVector3(1.0,0.0,-1.0).normalize());

	continuity = 2;
	support_radius = 5;
	surflet_degree = 2;
	param_changed = true;

	base_level_offset=1.0;

	updateRBF();
}

//immediately interpolate a particular set of data
SurfletCompactRBF::SurfletCompactRBF(int cont,double support,int surflet_deg,float bloffset,const std::vector<gmVector3> &positions, 
							 const std::vector<double> &interp_vals,
							 const std::vector<gmVector3> &normals)
{
	continuity= cont;
	support_radius=support;
	surflet_degree=surflet_deg;
	param_changed=true;
	base_level_offset=bloffset;

	std::valarray<bool> temp;
	interpolate(positions,interp_vals,normals,temp,true);
}

/** Deallocate Surflets
  */
void SurfletCompactRBF::clearSurflets()
{
	for(unsigned int c=0;c<surflets.size();c++)
	{
		if(surflets[c]) delete surflets[c];
	}
	surflets.clear();
}


/** Returns the number of parameters of the Surflet Compact RBF (which is 3).
 */
unsigned int SurfletCompactRBF::qlen()
{
	return 3;
}

/**
 * Loads parameters into q.
 * @param q An array of doubles, of size qlen, into which to copy the parameter array of the surface.
 */
void SurfletCompactRBF::getq(double *q)
{
	q[0] = static_cast<double>(continuity);
	q[1] = support_radius;
	q[2] = static_cast<double>(surflet_degree);
}

/** Stores parameters into q.
 * @param q An array of doubles, of size qlen, containing the values to store in the parameter vector.
 */
void SurfletCompactRBF::_setq(double *q)
{
	if(static_cast<int>(q[0])!=continuity) param_changed=true;
	if(support_radius!=q[1]) param_changed=true;
	if(static_cast<int>(q[2])!=surflet_degree) param_changed=true;
	continuity = static_cast<int>(q[0]);
	support_radius = q[1];
	surflet_degree = static_cast<int>(q[2]);
	
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
void SurfletCompactRBF::getqname(char **qn)
{
	qn[0] = "Continuity";
	qn[1] = "Support Radius";
	qn[2] = "Surflet Degree";
}


#ifndef INTERVAL_EVAL_ONLY
/** Evaluate the variational implicit surface.
 * The surface consists of a list of radial basis functions and a bunch of surflets 
 */
double SurfletCompactRBF::proc(const gmVector3 & x)
{
  double fVal = 0.0;
  
  std::vector<int> nz_centers;
  support_tree.radiusQuery(x, support_radius,nz_centers);
  
  // If the point is outside of the support radius of any point, make its value huge.
  //  This will cause ParticleValueDeath (if it is present) to destroy any floater
  //  outside of the support region, eliminating the "runaway floaters" problem.
/*  if(nz_centers.size() == 0)
	  return 100.0;*/

  for(unsigned int i = 0; i < nz_centers.size(); i++)
  {
	  int ci=nz_centers[i];
	  if(surflets[ci])
	  {
		  fVal += (surflets[ci]->eval(x) + centers[ci].w) * phi(x - centers[ci].c);
	  }else
	  {
		  fVal += (centers[ci].w) * phi(x - centers[ci].c);
	  }
  }

  return fVal+base_level_offset;
}

gmVector3 SurfletCompactRBF::grad(const gmVector3 & x)
{
//  return Implicit::grad(x);
  gmVector3 g(0.0,0.0,0.0);

  std::vector<int> nz_centers;
  support_tree.radiusQuery(x, support_radius,nz_centers);

  for (unsigned int i = 0; i < nz_centers.size(); i++)
  {
	  int ci=nz_centers[i];
	  if(surflets[ci])
	  {
		  //product rule!
			g += surflets[ci]->eval(x) * gradphi(x - centers[ci].c);
			g += surflets[ci]->gradient(x) * phi(x - centers[ci].c);
	  }

   	  g += centers[ci].w * gradphi(x - centers[ci].c);	  
  }

  return g;
}

inline gmMatrix3 tensorProduct(const gmVector3 &v1,const gmVector3 &v2)
{	
	return gmMatrix3(v1[0]*v2[0],v1[0]*v2[1],v1[0]*v2[2],
		             v1[1]*v2[0],v1[1]*v2[1],v1[1]*v2[2],
					 v1[2]*v2[0],v1[2]*v2[1],v1[2]*v2[2]);
}

gmMatrix3 SurfletCompactRBF::hess(const gmVector3 & x)
{
	return Implicit::hess(x);
  /*gmMatrix3 g;

  std::vector<int> nz_centers;
  support_tree.radiusQuery(x, support_radius,nz_centers);

  for (unsigned int i = 0; i < nz_centers.size(); i++)
  {
	  int ci=nz_centers[i];
	  if(surflets[ci])
	  {
		  //product rule!
		  g += surflets[ci]->eval(x) * hessphi(x - centers[ci].c) + 2.0*tensorProduct(surflets[ci]->gradient(x),gradphi(x - centers[ci].c));
		  g += surflets[ci]->hessian(x) * phi(x - centers[ci].c);

	  }
	  g += centers[ci].w * hessphi(x - centers[ci].c);	  
  }

  return g;*/
}

#endif

/** Update the surface
  */
void SurfletCompactRBF::updateRBF(void)
{
  int i,j,n;

  //converts constraints, builds KDtree, and fits surflets
  convert_constraints();
  n = centers.size();

  TNT::Matrix<double> A(n,n);
  TNT::Vector<int> ipiv(n);
  TNT::Vector<double> h(n);

  for(i=0; i < n; i++)
  {
	for(j=0; j < n; j++) 
	{
      A[i][j] = phi(centers[j].c - centers[i].c);
	}
  }

  for(i = 0; i < n; i++)
    h[i] = centers[i].h - base_level_offset;

  //Account for surflet contributions
  for(i=0;i<n;i++)
  {    
    for(j=0;j<n;j++)
    {
		if(surflets[j]) h[i]-=surflets[j]->eval(centers[i].c)*(A[i][j]);
    }
  }


  TNT::LU_factor(A,ipiv);
  TNT::LU_solve(A,ipiv,h);

  for(i = 0; i < n; i++)
    centers[i].w = h[i];
}

/** Convert the "conceptual" constraints into actual Surflet RBF constraint points (fit surflets)
 */      
void SurfletCompactRBF::convert_constraints()
{
	centers.clear();
	clearSurflets();

	int n = constraints.size();
	  // Insert all of the centers into the kd-tree.
	std::vector<gmVector3> points(n);
	for(int k = 0; k < n; k++)
	{
		points[k]=constraints[k].c;
	}
	support_tree.generate(points);

	// Convert all the constraint points.  Normal constraints get surflets fitted to them
	//  offset function value constraints.
	for(unsigned int i = 0; i < constraints.size(); i++)
	{
		centers.push_back(RBFPoint(constraints[i].c, constraints[i].h));
		if(constraints[i].nc)
		{
			//find all points near us
			std::vector<int> nz_centers;
			support_tree.radiusQuery(constraints[i].c, support_radius,nz_centers);

			surflets.push_back(newSurflet());
			//get weights
			std::vector<double> nz_weights(nz_centers.size(),1.0);
			for(unsigned int c=0;c<nz_centers.size();c++)
			{
				nz_weights[c]=phi(constraints[nz_centers[c]].c-constraints[i].c);
			}
			surflets.back()->fit(constraints[i].c,constraints[i].n,constraints,nz_centers,nz_weights);
		}else
		{
			surflets.push_back(0); //points that aren't normal constraints get no surflet
		}
	}
}

/** 
 * Loads control particles into RBF vector and calls updateRBF().
 * @param positions A vector of constraint positions.
 * @param interp_vals A vector of desired function values for each constraint (determines interior, exterior, or surface constraint).
 * @param normals A vector of the desired normals for normal constraints.
 * @param flexible A vector which serves no discernable purpose for RBFs.
 * @param pos_changed Have any of the constraints changed since we last solved for a surface?  If not, do not re-solve.
 */
void SurfletCompactRBF::interpolate(std::vector<gmVector3> positions, std::vector<double> interp_vals, std::vector<gmVector3> normals, std::valarray<bool>& flexible, bool pos_changed)
{
	if(pos_changed || param_changed)
	{
		constraints.clear();
		constraints.resize(positions.size());

		for(unsigned int i = 0; i < positions.size(); i++)
		{
			// Normal constraint
			if(normals[i] != gmVector3(0.0, 0.0, 0.0))
				constraints[i] = RBFModelerConstraint(positions[i], 
												normals[i],
												interp_vals[i], true);
			// Non-normal constraint
			else
				constraints[i] = RBFModelerConstraint(positions[i], normals[i], interp_vals[i], false);
		}
	
		updateRBF();
		param_changed = false;
	}
}

/*************************LINEARSURFLET IMPLEMENTATION*************/
void LinearSurflet::fit(gmVector3 point,gmVector3 normal,const RBFModelerConstraints &points,const std::vector<int> &nz_indices,const std::vector<double> &nz_weights)
{
	//for a linear surflet, there's not much fitting to be done, just find the plane equation that matches the normal and the point
	n=normal;
	n.normalize();
	d = (float) -dot(normal,point);
}

float LinearSurflet::eval(const gmVector3 &point) const
{
	return (float) dot(n,point)+d;
}

gmVector3 LinearSurflet::gradient(const gmVector3 &point) const
{
	return n; //gradient of a linear surflet is trivial
}

/*********************QUADRATICSURFLET IMPLEMENTATION*************/
//utility for choosing an orthogonal 3D basis given just one of the axes
void chooseOrthogonalBasis(gmVector3 &u,gmVector3 &v,const gmVector3 &knownZ)
{
	//It doesn't matter what two vectors we pick for the plane basis vectors
	//as long as they are perpendicular to the ray.
	//So, we pick the axis direction which the ray direction is smallest in and use it to produce
	//a vector that is truly perpendicular to knownZ
	int mini=0;		
	for(int i=1;i<3;i++) if(fabs(knownZ[i])<fabs(knownZ[mini])) mini=i;
	u=gmVector3(0.0,0.0,0.0);
	u[mini]=1.0;
	//subtract the portion that's parallel to the ray
	u-=knownZ*(dot(knownZ,u));
	//renormalize
	u.normalize();
	//produce the other basis vector (perp to both u and ray.direction)
	v=cross(knownZ,u);
}

void QuadraticSurflet::fit(gmVector3 point,gmVector3 normal,const RBFModelerConstraints &points,const std::vector<int> &nz_indices,const std::vector<double> &nz_weights)
{
	//produce an orthonormal coordinate system, with one axis along normal direction,
	//centered at point
	what=normal;
	chooseOrthogonalBasis(uhat,vhat,what);
	origin=point;

	//find out how many points are within the non-zero radius of support
	int n=nz_indices.size();
	if(n<12) //if there are less than twelve, revert to LinearSurflet behavior
	{
		//even though we technically need only six points to do this fit, if we have a small
		//number of points, the chances are higher that the matrix will be nearly singular,
		//causing various problems. Thus, we should not use quadratic fitting when the number of
		//points is very small.
		a=b=c=d=e=f=0.0;
		return;
	}

	TNT::Matrix<double> A(n,6);
	TNT::Vector<double> rhs(n);

	//build rectangular matrix
	for(int j=0;j<n;j++)
	{
		gmVector3 local=toLocal(points[nz_indices[j]].c);
		double weight = sqrt(nz_weights[j]); 
		A[j][0]=local[0]*local[0]*weight;
		A[j][1]=2.0*local[0]*local[1]*weight;
		A[j][2]=local[1]*local[1]*weight;
		A[j][3]=local[0]*weight;
		A[j][4]=local[1]*weight;
		A[j][5]=weight;
		rhs[j]=local[2]*weight;
	}

	//perform QR factorization
	JAMA::QR<double> QRfactors(A);
	if(QRfactors.isFullRank())
	{
		TNT::Vector<double> x(6);
		x=QRfactors.solve(rhs);
		//the solve succeeded! gather parameters
		a=(float)x[0];
		b=(float)x[1];
		c=(float)x[2];
		d=(float)x[3];
		e=(float)x[4];
		f=(float)x[5];
	}else
	{
		//matrix was singular, revert to LinearSurflet behavior
		a=b=c=d=e=f=0.0;
	}
}

float QuadraticSurflet::eval(const gmVector3 &point) const
{
	gmVector3 local=toLocal(point);
	return (float)(local[2] - (a*local[0]*local[0] + 2.0*b*local[0]*local[1] + c*local[1]*local[1] + d*local[0] + e*local[1] + f));
}

gmVector3 QuadraticSurflet::gradient(const gmVector3 &point) const
{
	gmVector3 local=toLocal(point);
	//compute gradient in local frame
	gmVector3 local_g(-2.0*a*local[0]-2.0*b*local[1]-d,
					  -2.0*b*local[0]-2.0*c*local[1]-e,
					  1.0);
	//convert to global frame and return
	return local_g[0]*uhat+local_g[1]*vhat+local_g[2]*what;
}

gmMatrix3 QuadraticSurflet::hessian(const gmVector3 &point) const
{
	//THIS DOES NOT WORK YET!
	//But that is OK, because I am currently using finite differences for the Hessian
	return gmMatrix3();
}
