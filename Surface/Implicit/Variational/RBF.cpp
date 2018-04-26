/** @file RBF.cpp 
 * Implementation of variational implicit surfaces.
 * @author William Nagel
 * @date Fall Semester, 2000.
 */
//procq was overloaded to return an array of 0s... thus
//we avoid getting setq called all the time with different parameters to
//do a meaningless finite-differences computations in child classes. -Scott Kircher

#include "RBF.h"

#include "Particles/Particles.h"
#include "Particles/Attributes/ParticlePosition.h"
#include "Particles/Attributes/ParticleScalar.h"
#include "Particles/Attributes/ParticleDensity.h"

#include "tnt/tnt.h"
#include "tnt/cmat.h"
#include "tnt/vec.h"
#include "tnt/lu.h"
#include <stdio.h>


// implicit should not depend on particles, the surface should have its points stored in its own format
// particles are a way of displaying surfaces. This is a design flaw.
//#include "ParticlesInclude.h"
#include "Surface/gmTNTconvert.h"
#include "Particles/Attributes/ParticlePosition.h"
#include "Particles/Attributes/ParticleScalar.h"

REGISTER_IMPLICIT(RBF,"Variational:RBF");

/** Default empty constructor. This currently builds an RBF sphere consisting of six 0 centers
 * at (+/-1, +/-1, +/-1) and one -1 center at (0,0,0).
 */


RBF::RBF(void)
{
#if 0
	p = NULL;
	new SurfParRefParam(this,&p,"<empty>","pCenters","Centers",
		"Collection of particle whose positions provide the center locations for the radial basis functions.");
#else
	new SurfAttrRefParam(this,(ParticleAttribute **)&positions,"ParticlePosition","pos","centers",
		"Particle position attribute used for RBF center location.");
	new SurfAttrRefParam(this,(ParticleAttribute **)&values,"values","val","values",
		"Desired scalar value at the particle position.");
// this parameter is unused...
//	new SurfAttrRefParam(this,(ParticleAttribute **)&weights,"weights","weights","weights",
//		"Weights of RBF functions cached per-particle.");
#endif
	//new SurfStringParam();
	//new SurfStringParam(Surface *parent, std::string s, std::string sn, std::string n, std::string d);
	//new SurfStringParam(this, filename, "rbfFilename","rbf value filename","file from which the rbf center values are read.");
	//new SurfStringParam(this,&filename,"<empty>","rbfFile","rbf Values File","file from which to load the values of the rbf constraint points");


	constraints.resize(7);

	/*constraints[0] = RBFModelerConstraint(gmVector3(1.0,0.0,0.0), gmVector3(1.0,0.0,0.0));
	constraints[1] = RBFModelerConstraint(gmVector3(0.0,1.0,0.0), gmVector3(0.0,1.0,0.0));
	constraints[2] = RBFModelerConstraint(gmVector3(0.0,0.0,1.0), gmVector3(0.0,0.0,1.0));
	constraints[3] = RBFModelerConstraint(gmVector3(-1.0,0.0,0.0), gmVector3(-1.0,0.0,0.0));
	constraints[4] = RBFModelerConstraint(gmVector3(0.0,-1.0,0.0), gmVector3(0.0,-1.0,0.0));
	constraints[5] = RBFModelerConstraint(gmVector3(0.0,0.0,-1.0), gmVector3(0.0,0.0,-1.0));*/

	/*constraints[0] = RBFModelerConstraint(gmVector3(1.0,0.0,0.0),0);
	constraints[1] = RBFModelerConstraint(gmVector3(0.0,1.0,0.0),0);
	constraints[2] = RBFModelerConstraint(gmVector3(0.0,0.0,1.0),0);
	constraints[3] = RBFModelerConstraint(gmVector3(-1.0,0.0,0.0),0);
	constraints[4] = RBFModelerConstraint(gmVector3(0.0,-1.0,0.0),0);
	constraints[5] = RBFModelerConstraint(gmVector3(0.0,0.0,-1.0),0);
	constraints[6] = RBFModelerConstraint(gmVector3(0.0,0.0,0.0),-2);*/

	centers.resize(7);
	centers[0] = RBFPoint(gmVector3(1.0,0.0,0.0),0);
	centers[1] = RBFPoint(gmVector3(-1.0,0.0,0.0),0);
	centers[2] = RBFPoint(gmVector3(0.0,1.0,0.0),0);
	centers[3] = RBFPoint(gmVector3(0.0,-1.0,0.0),0);
	centers[4] = RBFPoint(gmVector3(0.0,0.0,1.0),0);
	centers[5] = RBFPoint(gmVector3(0.0,0.0,-1.0),0);
	centers[6] = RBFPoint(gmVector3(0.0,0.0,0.0),-1);

	param_changed = true;
	updateRBF();
}


//[TK 5.18.05]
//SurfParamTextCtrl::OnValue is called when a name is assigned to a child of binaryOp
// then it calls surfaces.attach()


// overridden to attach a particle system to use as the centers for this RBF object
// this is called from surfaceTree.onvalue(...) where parameters are updated when new values are entered through the GUI

void RBF::_setq(double* q)
{
	int i, n = qlen();

//	centers.resize(n/3+1);
	centers.resize(n/3);
	for (i=0; i < n/3; i++)
	{
		centers[i].c[0] = q[i*3];
		centers[i].c[1] = q[i*3+1];
		centers[i].c[2] = q[i*3+2];
//		centers[i].w = q[i*4+3];
//		centers[i].w = q[i];
		if ( values )
			centers[i].h = values->getScalar(i);
		if ( positions )	
			positions->setPosition(i, centers[i].c);
//		centers[i].h = 0;
	}
	//for (int i=0; i<m_numCoef; i++) 
    //m_a[i] = q[i];
	// attach a particle system.
	// load the verticies from that particle system as the centers of this implicit


	/*
	m_p[0] = q[i*4];
	m_p[1] = q[i*4+1];
	m_p[2] = q[i*4+2];
	m_p[3] = q[i*4+3];
	*/

//	centers[i].c[0] = 0;
//	centers[i].c[1] = 0;
//	centers[i].c[2] = 0;
//	centers[i].h = -1;

	updateRBF();
}

unsigned int RBF::qlen()
{
	return centers.size() * 3; // + 4;
}

void RBF::getq(double *q)
{
	int i, n = qlen(); // - 4;

	for (i=0; i < n/3; i++)
	{
		q[i*3] = centers[i].c[0];
		q[i*3+1] = centers[i].c[1];
		q[i*3+2] = centers[i].c[2];
//		q[i*4+3] = centers[i].w;
//		q[i] = centers[i].w;
	}

	/*q[i*4] = m_p[0];
	q[i*4+1] = m_p[1];
	q[i*4+2] = m_p[2];
	q[i*4+3] = m_p[3];*/
}


void RBF::getqname(char **qn)
{
	char coord[] = { 'x', 'y', 'z' };
	for (unsigned int i=0; i < qlen(); i++) {
		qn[i] = new char[15];
		sprintf(qn[i], "center %d.%c", (i / 3)+1, coord[i%3]);
	}
//	qn[0] = "pCenters";
}


//overridden to load constraints from the file
bool RBF::readImplicit(std::ifstream &file,bool verbose)
{
	std::vector<double> params;
	readParameters(file,params);

	int num_values, starting_point;
	int numq=qlen();
	// If this is an RBF, the rest of the values represent control point locations, normals, and
	//  function values (0 = surface, -1 = interior, +1 = exterior).

	// Some of our child classes have actual parameters, so check qlen to see how many
	// we have
	// Put in the real parameters.
	if(numq>0)
	{
		std::vector<double> real_params = std::vector<double>(numq, 0.0);
		for(unsigned int c=0;c<(unsigned int) (numq) && c<params.size();c++)
			real_params[c] = params[c];

		setq(real_params);
			
		num_values = params.size() - numq;
		starting_point = numq;
	}
	else
	{
		num_values = params.size();
		starting_point = 0;
	}

	// Pull out the rest of the parameters 7 at a time and create constraints out of them.
	//  If the number of values is not a multiple of 7, throw away whatever is left.
	if(num_values % 7)
		num_values -= (num_values % 7);

	unsigned int j = starting_point;
	constraints.clear();

	while(j < params.size())
	{
		gmVector3 position, normal;
		double fun_val;
		bool normal_constraint;

		position = gmVector3(params[j], params[j + 1], params[j + 2]);
		normal = gmVector3(params[j + 3], params[j + 4], params[j + 5]);
		fun_val = params[j + 6];
		
		if(normal.lengthSquared() == 0.0)
			normal_constraint = false;
		else
			normal_constraint = true;

		constraints.push_back(RBFModelerConstraint(position, normal, fun_val, normal_constraint));

		j += 7;
	}
	updateRBF();

	return true;
}

//void RBF::interpolate(Particles *ps, float phi)
void RBF::interpolate()
{
	if (!positions || !values)
		return;
	centers.resize(positions->x.size());

	for(unsigned int i = 0; i < positions->x.size(); i++) {
//		double a = positions->x[i][0];
//		double b = positions->x[i][1];
//		double c = positions->x[i][2];
//		std::cout << "Center " << i << ": " << positions->x[i][0] << "," << positions->x[i][1] << "," << positions->x[i][2] << std::endl;
		centers[i] = RBFPoint(positions->x[i], values->getScalar(i));
	}

	// add a negative center..
//	centers[positions->x.size()] = RBFPoint(gmVector3(0, 0, 0), -1);
	
	updateRBF();
}

void RBF::updateRBF(void)
{
  int i,j,n;

  n = centers.size();

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

  int singularTest = TNT::LU_factor(A,ipiv); // if the matrix is singular, this call will break the program [TK 3.05]
  if(singularTest != 0)
	  return; // the system is singular, and cannot be solved.  That would crash the program immediately.
  TNT::LU_solve(A,ipiv,h);

  for(i = 0; i < n; i++)
    centers[i].w = h[i];

  m_p[0] = h[n];
  m_p[1] = h[n+1];
  m_p[2] = h[n+2];
  m_p[3] = h[n+3];
  param_changed = false;
}

#ifndef INTERVAL_EVAL_ONLY
/** Evaluate the variational implicit surface.
 * The surface consists of a list of radial basis functions and
 * a plane (the affine part).
 */

void RBF::procq(const gmVector3 & x, double *dq)
{
	int i, n = centers.size();

	// dq should be size 3*n (each x,y,z component gets an entry)
	for (i=0; i < n; i++)
	{
		gmVector3 r = x - centers[i].c;

		double coeff = -3.0 * r.length() * centers[i].w;
		dq[i*3] = coeff * (x[0] - centers[i].c[0]);
		dq[i*3+1] = coeff * (x[1] - centers[i].c[1]);
		dq[i*3+2] = coeff * (x[2] - centers[i].c[2]);

//		dq[i*4+3] = phi(x - centers[i].c);

//		dq[i] = phi(x - centers[i].c);
	}

//	dq[i*4] = 1;
//	dq[i*4+1] = x[0];
//	dq[i*4+2] = x[1];
//	dq[i*4+3] = x[2];
}


double RBF::proc(const gmVector3 & x)
{
  double fVal = 0.0;
  int numC = (int)centers.size();
  for(unsigned int i = 0; i < (unsigned int)numC; i++)
    fVal += centers[i].w * phi(x - centers[i].c);

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
gmVector3 RBF::grad(const gmVector3 & x)
{
  gmVector3 g(0.0,0.0,0.0);

  for (unsigned int i = 0; i < centers.size(); i++)
	g += centers[i].w * gradphi(x - centers[i].c);

  g += gmVector3(m_p[1],m_p[2],m_p[3]);

  return g;
}

gmMatrix3 RBF::hess(const gmVector3 & x)
{
  gmMatrix3 g;

  for (unsigned int i = 0; i < centers.size(); i++)
    g += centers[i].w * hessphi(x - centers[i].c);

  return g;
}
#endif

/** 
 * Evaluate the variational implicit surface.  This is the
 * Interval version of proc.
 * @param b The Box3d region over which to evaluate the function.
 */
Intervald RBF::proc(const Box<double>&  b) 
{
  Intervald fVal(0.0);
  
  for (unsigned int i = 0; i < centers.size(); i++)
    fVal += Intervald(centers[i].w) * phi(Box3d(b - Box3d(centers[i].c)));

  fVal += Intervald(m_p[3]) * b[2];
  fVal += Intervald(m_p[2]) * b[1];
  fVal += Intervald(m_p[1]) * b[0];
  fVal += Intervald(m_p[0]);

  return fVal; 
};

/**
 * Gradient of the sum of radial basis functions.
 * This is the Intefvald/Box3d version.
 */
Box3d RBF::grad(const Box<double>&  b)
{
  Box3d g;

  for (unsigned int i = 0; i < centers.size(); i++)
    g = g + gradphi(Box3d(b - Box3d(centers[i].c))) * Intervald(centers[i].w);

  g = g + Box3d(Intervald(m_p[1]),Intervald(m_p[2]),Intervald(m_p[3]));

  return g;
}

IMatrix3d RBF::hess(const Box<double>&  b)
{
  IMatrix3d g;
  IMatrix3d temp;

  for (unsigned int i = 0; i < centers.size(); i++)
    g = g + hessphi(Box3d(b - Box3d(centers[i].c))) * Intervald(centers[i].w);

  return g;
}

/** Evaluate the radial basis function.
 * For 3-D, the radial basis function is ||x||^3.
 */
double RBF::phi(const gmVector3 & x)
{
  return pow(x.lengthSquared(), 1.5);
}

/** Evaluate the radial basis function over an interval.
 * @param x The interval over with to evaluate ||x||^3
 */
Intervald RBF::phi(const Box<double>&  b)
{
  Box3d b1 = Box3d(b[0],b[1],b[2]);
  return b1.length().pow(3);
}

gmVector3 RBF::gradphi(const gmVector3 & x)
{
  double r = x.length();
  return 3.0*r*x;
}

Box3d RBF::gradphi(const Box<double>&  b)
{
  Box3d b1 = Box3d(b[0],b[1],b[2]);
  b1 *= b.length();
  b1 *= Intervald(3.0);
  return b1;
}

gmMatrix3 RBF::hessphi(const gmVector3 & x)
{
  double length = x.length();

  // Cheat - added 6/5/04 by Mike Flavin
  //  - If a silhouette particle is in the same location as a control particle
  //    (which it always will be, since all particles in RBFs start at the control points),
  //    the program will crash because of divide by zero here.  So, I'm cheating by making
  //    sure length is NEVER zero.  Hopefully this will work.
  if(gmIsZero(length))
	  length = 0.000001;

  return 3.0 * 
    (gmMatrix3::identity() * x.lengthSquared() + outer(x,x)) /
      length;
}

IMatrix3d RBF::hessphi(const Box<double>&  b)
{ 
  int i,j;
  IMatrix3d bsq;
  IMatrix3d sum;
  Box3d b1 = Box3d(b[0],b[1],b[2]);

  bsq[0][0] = bsq[1][1] = bsq[2][2] = b1.lengthSquared();
  sum = IMatrix(bsq) + IMatrix(outer(b1,b1));
  for (i = 0; i < sum.num_rows(); i++)
    for (j = 0; j < sum.num_cols(); j++)
      sum[i][j] *= (Intervald(3.0) / b1.length());

  return sum;
}

