/**
 * Implementation of a rbf thin plate spline (r^2 log r)
 * @file ThinPlateSpline.cpp
 * @date 10 April, 2006
 * @author Eric Lorimer
 * @remarks
 */

#include "ThinPlateSpline.h"
#include "Particles/Attributes/ParticlePosition.h"
#include "Particles/Attributes/ParticleScalar.h"

#include "tnt/tnt.h"
#include "tnt/cmat.h"
#include "tnt/vec.h"
#include "tnt/lu.h"


REGISTER_IMPLICIT(ThinPlateSpline,"ThinPlateSpline");

ThinPlateSpline::ThinPlateSpline(void)
{
	new SurfAttrRefParam(this,(ParticleAttribute **)&positions,"ParticlePosition","pos","centers",
		"Particle position attribute used for RBF center location.");
	new SurfAttrRefParam(this,(ParticleAttribute **)&values,"values","val","values",
		"Desired scalar value at the particle position.");

	centers.push_back( gmVector3(0, 0, 0) );
	centers.push_back( gmVector3(0, 1, 0) );
	centers.push_back( gmVector3(2, 0, 0) );

	weights.push_back( 0.9 );
	weights.push_back( -1.3 );
	weights.push_back( -0.3 );
}


double ThinPlateSpline::proc(const gmVector3 & x)
{
	double minp = 1e20;
	for (unsigned int i=0; i < centers.size(); i++)
		if ( proci(x, i) < minp )
			minp = proci(x, i);

	return minp;
}


double ThinPlateSpline::proci(const gmVector3 & x, int i)
{
	gmVector3 center = centers[i];
	double xd = x[0] - center[0];
	double yd = x[1] - center[1];
	double r = sqrt(xd * xd + yd * yd);
	if ( weights[i] > 0)
		return weights[i] * r*r*log(r) - x[2];
//	if ( weights[i] <= 0)
	return x[2] - weights[i] * r*r*log(r);
}


gmVector3 ThinPlateSpline::gradi(const gmVector3 & x, int i)
{
	gmVector3 gradient;
	gmVector3 center = centers[i];

	double xd = x[0] - center[0];
	double yd = x[1] - center[1];

	gradient[0] = weights[i] * xd * (log(xd*xd + yd*yd) + 1);
	gradient[1] = weights[i] * yd * (log(xd*xd + yd*yd) + 1);
	gradient[2] = -1;

	if ( weights[i] > 0 ) return gradient;
	if ( weights[i] < 0 ) return -gradient;
	return gradient;
}


gmVector3 ThinPlateSpline::grad(const gmVector3 & x)
{
	unsigned int i, mi = 0;
	double px = 1e20;
	for (i=0; i < centers.size(); i++) {
		double f = proci(x, i);
		mi = (f < px) ? i : mi;
		px = (f < px) ? f : px;
	}

	return gradi(x, mi);
}


void ThinPlateSpline::interpolate(Particles *ps, float phi)
{
	if (!positions || !values)
		return;
	centers.resize(positions->x.size());

	for(unsigned int i = 0; i < positions->x.size(); i++)
		centers[i] = positions->x[i];

	updateRBF();
}



void ThinPlateSpline::updateRBF(void)
{
  int i,j,n;

  n = centers.size();

  TNT::Matrix<double> A(n + 3, n + 3);
  TNT::Vector<int> ipiv(n + 3);
  TNT::Vector<double> h(n + 3);

  for(i=0; i < n; i++) {
	for(j=0; j < n; j++)
      A[i][j] = phi(centers[j] - centers[i]);

    A[n][i] = 1.0;
    A[i][n] = 1.0;

    for(j=0; j < 2; j++)
    {
      A[n + 1 + j][i] = centers[i][j];
      A[i][n + 1 + j] = centers[i][j];
    }
  }

  for(i = 0; i < n; i++)
	h[i] = values->getScalar(i);

  h[n] = 0.0;
  h[n + 1] = 0.0;
  h[n + 2] = 0.0;

  int singularTest = TNT::LU_factor(A,ipiv); // if the matrix is singular, this call will break the program [TK 3.05]
  if(singularTest != 0)
	  return; // the system is singular, and cannot be solved.  That would crash the program immediately.
  TNT::LU_solve(A,ipiv,h);

  for(i = 0; i < n; i++)
	weights[i] = h[i];

  m_p[0] = h[n];
  m_p[1] = h[n+1];
  m_p[2] = h[n+2];
}


void ThinPlateSpline::getq(double *q)
{
	for (unsigned int i=0; i < centers.size(); i++) {
		q[i*4] = weights[i];
		q[i*4+1] = centers[i][0];
		q[i*4+2] = centers[i][1];
		q[i*4+3] = centers[i][2];
	}
}


void ThinPlateSpline::_setq(double *q)
{
	for (unsigned int i=0; i < centers.size(); i++) {
		weights[i] = q[i*4];
		centers[i][0] = q[i*4+1];
		centers[i][1] = q[i*4+2];
		centers[i][2] = q[i*4+3];
	}
}


void ThinPlateSpline::getqname(char **qn)
{
	for (unsigned int i=0; i < centers.size(); i++) {
		qn[i*4] = "Weight";
		qn[i*4+1] = "Center X";
		qn[i*4+2] = "Center Y";
		qn[i*4+3] = "Center Z";
	}
}


double ThinPlateSpline::phi(const gmVector3 & x)
{
	double r = sqrt(x[0] * x[0] + x[1] * x[1]);
	if ( r == 0 ) return 0;
	return r*r * log(r);
}