/** @file RBF.h
 * Declaration of an implicit surface class based on radial basis function interpolation.
 * It offers the possibility to exchange the function phi (important for smoothing Carr et al. 2003). 
 * Also calculation of the RBF can done based on particles
 * Therefore it is possible to create a particles object whose position and values can be connected.
 * E.g.: Create a particle system that spreads out over the surface based on curvature and then use these 
 * particles to calculate an appropriate RBF. Or place particles on a mesh etc.
 * @author: Elmar Eisemann and I stole a lot ! ;)
 * @date: Sommer, 2006
**/
#include "InterfaceRBF.h"
#include "RBFBasicFunctionsAndSelector.h"

#include "Particles/Particles.h"
#include "Particles/Attributes/ParticlePosition.h"
#include "Particles/Attributes/ParticleScalar.h"
#include "Particles/Attributes/ParticleDensity.h"

#include "tnt/tnt.h"
#include "tnt/cmat.h"
#include "tnt/vec.h"
#include "tnt/lu.h"
#include <stdio.h>

#include "Surface/gmTNTconvert.h"
#include "Particles/Attributes/ParticlePosition.h"
#include "Particles/Attributes/ParticleScalar.h"

REGISTER_IMPLICIT(InterfaceRBF,"Variational:InterfaceRBF");

/** Default empty constructor. This currently builds an RBF sphere consisting of six 0 centers
 * at (+/-1, +/-1, +/-1) and one -1 center at (0,0,0).
 */


InterfaceRBF::InterfaceRBF(void)
:RBF() 
{
	new SurfAttrRefParam(this,(ParticleAttribute **)&zeroConstraints,"ZeroConstraints","zeroPos","ZeroCons",
		"Particle positions used as zero constraints.");

	new SurfParamButton(this, new InterpolateParticlesWithRBF(this),"interpolate particles","fit RBF",
		"find RBF corresponding to particle positions and values");
	//I do not really like this kind of construction
	//but as it is done everywhere I stick to it...
	new	SurfParamComboBox(this,new RBFBasicFunctionSelector(this, &_phiFunction), "RBFSelect","RBFSelector",
		"The selection decides which basic functions to use during creation and evaluation. Note: The creation choice does not fix the evaluation choice.");
	new SurfParamString(this,&_filename,"","fn","filename","Load/Save filename (including location)");
	new SurfParamButton(this, new LoadRBFCallback(this),"loadRBF", "load solution", "loadRBF from file");
	new SurfParamButton(this, new SaveRBFCallback(this), "saveRBF", "save solution", "saveRBF to file");


	//Init with default RBF
	centers.resize(7);
	centers[0] = RBFPoint(gmVector3(1.0,0.0,0.0),0);
	centers[1] = RBFPoint(gmVector3(-1.0,0.0,0.0),0);
	centers[2] = RBFPoint(gmVector3(0.0,1.0,0.0),0);
	centers[3] = RBFPoint(gmVector3(0.0,-1.0,0.0),0);
	centers[4] = RBFPoint(gmVector3(0.0,0.0,1.0),0);
	centers[5] = RBFPoint(gmVector3(0.0,0.0,-1.0),0);
	centers[6] = RBFPoint(gmVector3(0.0,0.0,0.0),-1);
	updateRBF();
}

void InterfaceRBF::load()
{
	std::ifstream file;
	file.open(_filename.c_str(),std::ios::in | std::ios::binary);
	if (file.is_open())
	{
		//load RBF
		unsigned nbCenters;
		file.read(reinterpret_cast<char*>(&nbCenters),sizeof(unsigned int));
		file.read(reinterpret_cast<char*>(m_p),sizeof(double)*4);

		centers.resize(nbCenters);
		for (unsigned int i=0;i<nbCenters;++i)
		{
			file.read(reinterpret_cast<char*>(&(centers[i].c[0])),sizeof(double));
			file.read(reinterpret_cast<char*>(&(centers[i].c[1])),sizeof(double));
			file.read(reinterpret_cast<char*>(&(centers[i].c[2])),sizeof(double));
			file.read(reinterpret_cast<char*>(&(centers[i].w)),sizeof(double));
		}
		file.close();
	}
	else
	{
		std::cout<<"InterfaceRBF:FileNotFound!\n";
	}
}
void InterfaceRBF::save() const
{
	std::ofstream file;
	file.open (_filename.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
	
	if (file.is_open())
	{
		//load RBF
		unsigned int size=centers.size();
		file.write(reinterpret_cast<const char*>(&size),sizeof(unsigned int));
		file.write(reinterpret_cast<const char*>(m_p),sizeof(double)*4);
		for (unsigned int i=0;i<centers.size();++i)
		{
			file.write(reinterpret_cast<const char*>(&(centers[i].c[0])),sizeof(double));
			file.write(reinterpret_cast<const char*>(&(centers[i].c[1])),sizeof(double));
			file.write(reinterpret_cast<const char*>(&(centers[i].c[2])),sizeof(double));
			file.write(reinterpret_cast<const char*>(&(centers[i].w)),sizeof(double));
		}
		file.close();
	}
	else
	{
		std::cout<<"Could not open file for writing!\n";
	}

}

void InterfaceRBF::_setq(double* q)
{
}

unsigned int InterfaceRBF::qlen()
{
	return 0;
}

void InterfaceRBF::getq(double *q)
{
}


void InterfaceRBF::getqname(char **qn)
{
}

void InterfaceRBF::procq(const gmVector3 & x,double *dq)
{

}



bool InterfaceRBF::readImplicit(std::ifstream &file,bool verbose)
{
	//this needs to be adjusted one day... 
	return RBF::readImplicit(file,verbose);
}




void InterfaceRBF::interpolate()
{
	if (!positions || !values)
		return;
	centers.resize(positions->x.size());

	for(unsigned int i = 0; i < positions->x.size(); i++) {
		centers[i] = RBFPoint(positions->x[i], values->getScalar(i));
	}

	if (dynamic_cast<ParticlePosition*>((ParticleAttribute*)zeroConstraints))
	{
		for(unsigned int i = 0; i < zeroConstraints->x.size(); i++) 
		{
			centers.push_back(RBFPoint(zeroConstraints->x[i], 0));	
		}
	}


	// add a negative center..	
	updateRBF();
}

void InterfaceRBF::updateRBF(void)
{
  int i,j,n;

  n = centers.size();

  TNT::Matrix<double> A(n + 4, n + 4);
  TNT::Vector<int> ipiv(n + 4);
  TNT::Vector<double> h(n + 4);

  for(i=0; i < n; i++) {
	for(j=0; j < n; j++) {
      // A[i][j] = (*this)[j].phi(convert((*this)[i].c));
      A[i][j] = _phiFunction->phi(centers[j].c - centers[i].c);
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
}




/** Evaluate the variational implicit surface.
 * The surface consists of a list of radial basis functions and
 * a plane (the affine part).
 */



double InterfaceRBF::proc(const gmVector3 & x)
{
  double fVal = 0.0;
  unsigned int numC = centers.size();
  for(unsigned int i = 0; i < numC; i++)
    fVal += centers[i].w * _phiFunction->phi(x - centers[i].c);

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
gmVector3 InterfaceRBF::grad(const gmVector3 & x)
{
  gmVector3 g(0.0,0.0,0.0);

  for (unsigned int i = 0; i < centers.size(); i++)
	g += centers[i].w * _phiFunction->gradphi(x - centers[i].c);

  g += gmVector3(m_p[1],m_p[2],m_p[3]);

  return g;
}

gmMatrix3 InterfaceRBF::hess(const gmVector3 & x)
{
  gmMatrix3 g;

  for (unsigned int i = 0; i < centers.size(); i++)
    g += centers[i].w * _phiFunction->hessphi(x - centers[i].c);

  return g;
}

/** 
 * Evaluate the variational implicit surface.  This is the
 * Interval version of proc.
 * @param b The Box3d region over which to evaluate the function.
 */
Intervald InterfaceRBF::proc(const Box<double>&  b) 
{
  Intervald fVal(0.0);
  
  for (unsigned int i = 0; i < centers.size(); i++)
    fVal += Intervald(centers[i].w) * _phiFunction->phi(Box3d(b - Box3d(centers[i].c)));

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
Box3d InterfaceRBF::grad(const Box<double>&  b)
{
  Box3d g;

  for (unsigned int i = 0; i < centers.size(); i++)
    g = g + _phiFunction->gradphi(Box3d(b - Box3d(centers[i].c))) * Intervald(centers[i].w);

  g = g + Box3d(Intervald(m_p[1]),Intervald(m_p[2]),Intervald(m_p[3]));

  return g;
}

IMatrix3d InterfaceRBF::hess(const Box<double>&  b)
{
  IMatrix3d g;
  IMatrix3d temp;

  for (unsigned int i = 0; i < centers.size(); i++)
	  g = g + _phiFunction->hessphi(Box3d(b - Box3d(centers[i].c))) * Intervald(centers[i].w);

  return g;
}

