/**
 * @file Implicit.cpp
 *
 * There is a compiler-time define you can set called INTERVAL_EVAL_ONLY.
 * When this define is NOT set, everything behaves as it usually does.  That
 * means that proc(), grad(), and hess() need to be defined for both
 * gmVector3 (non-Intervals) and Box<double> (Intervals), and the two sets
 * of functions are independent of each other.  This is probably the
 * behavior you want most often since the gmVector3 versions of proc(),
 * grad(), and hess() may be faster than the Box<double> versions.
 *
 * However, when you define INTERVAL_EVAL_ONLY, the gmVector3 versions of
 * proc(), grad(), and hess() depend on the Box<double> versions.  They
 * basically change the passed-in gmVector3 into a Box3d, call the Interval
 * version of proc(), grad(), or hess(), and then return the center() of the
 * result, changing the answer into the appropriate non-Interval type.
 * Thus, you only need to define Interval versions of proc(), grad(), and
 * hess().  This is useful in testing to verify that the Interval versions
 * of the functions are equivalent to the non-interval versions.
 *
 * There's also a new change for the Interval versions of proc(), grad(),
 * and hess().  It used to be that the default versions (as defined in
 * Implicit) simply returned a zero object (ie. Intervald(0) for proc(),
 * Box3d(0) for grad(), and IMatrix3d(0) for hess()).  This meant that to
 * use the Interval versions of proc(), grad(), and hess(), you HAD to
 * define these functions in all subclasses.  As of February 6, 2003, if you
 * do NOT have INTERVAL_EVAL_ONLY defined, the default Interval versions now
 * call the non-interval versions.  In effect, the default Interval versions
 * of proc(), grad(), and hess() return the same result as the non-Interval
 * versions, and in fact mostly ignore the fact that you called these
 * functions with Intervals (ie. use the center() of the Intervals).
 */

#include "Surface.h"
#include "Implicit.h"
#include "Variational/CompactRBF.h"

#include "Particles/Particles.h"
#include "Particles/Attributes/ParticlePosition.h"
#include "Particles/Attributes/ParticleVelocity.h"
#include "Particles/svd.h"
#include <iostream>

#include "Surface/gmTNTconvert.h"
#include "tnt/cholesky.h"
#include "tnt/transv.h"
#include "tnt/trisolve.h"



/// Default Constructor
Implicit::Implicit()
:Surface()
{
  m_epsilon = 0.0001;
  qold = NULL;
  qoldsize = 0;
  qsChangedFlag = false;
}

/// Destructor
Implicit::~Implicit()
{
  delete [] qold;
}

/** Reads implicit parameters from a file.
 * The default implementation grabs first the number of parameters
 * and then the individual parameters (as doubles). Child classes
 * may override this to provide any kind of detailed loading capabilities
 * necessary. This is called immediately after the class's name is
 * discovered in the file, and the class is instantiated. Implementations
 * of this function should only consume exactly the amount of input
 * necessary. Implicit hierarchy construction is still handled in
 * ImpFileManager::readImplicit()
 * \param file The input file stream.
 * \param verbose Whether or not the function should output error messages to stderr.
 * \return Whether or not the load succeeded.
 * \author Scott Kircher
 */
bool Implicit::readImplicit(std::ifstream &file,bool verbose)
{
	std::vector<double> params;
	readParameters(file,params);

	//by default, params must be same length as qlen.
	//Child classes may override this behavior to do anything they wish.
	if(params.size() != qlen())
	{
		params.resize(qlen(),0.0);
	}
	setq(params);

	return true;
}

/** Common param loading.
 * Many child classes may need to have initial behavior similar to Implicit::readImplicit()
 * and can use this utility function to read an array of doubles from the file
 * they may then do whatever they wish with it (readImplicit() just calls setq using this array).
 * This function reads first the number of parameters to load, and then the actual paramters.
 * \param file The input file stream.
 * \param params Vector of doubles that will be filled with the paramters read from the file.
 * \author Scott Kircher
 */
void Implicit::readParameters(std::ifstream &file,std::vector<double> &params)
{
	int num_params = 0;
	file >> num_params;
	params.resize(num_params,0.0);

	for (int i = 0; i < num_params; i++)
	{
		file >> params[i];
	}
}


/**
 * proc returns the evaluation of the function that defines the Implicit
 * surface.  In normal circumstances, the gmVector3 version of proc defaults
 * to a pure virtual function which means that you must define
 * proc(gmVector3) in its subclasses.  This interval version defaults to
 * returning a zero-width interval version of the double proc(), in case
 * an interval version has not been defined.
 * @param x Where to evaluate the function.
 * @return The function result.  The function is positive outside the
 *         surface and negative inside.
 */
Intervald Implicit::proc(const Box<double>& x) 
{ 
  return Intervald(proc(x.center()));
}

/// Shortcut evaluator/operator for the proc() of an Implicit
double Implicit::operator()(const gmVector3 & x) 
{ 
  return proc(x); 
}

/// Shortcut evaluator/operator for the proc() of an Implicit
Intervald Implicit::operator()(const Box<double>& x) 
{ 
  return proc(x); 
}

/** 
 * grad returns the gradient of the function (proc()) for a point or region.
 * Basically, the gradient returns a vector of partial derivatives, or more
 * formally grad(proc()) = (dproc/dx,dproc/dy,dproc/dz).
 * @note Unless grad is overwritten, it defaults to the numerical
 * approximation of the gradient = (proc(x+eps) - proc(x))/eps ...
 * @param x Where to evaluate the gradient of the function.
 * @return The gradient result.
 * @sa setEpsilon(double)
 */
gmVector3 Implicit::grad(const gmVector3 & x)
{
  gmVector3 gradient;
  gradient[0] = Fx(x);
  gradient[1] = Fy(x);
  gradient[2] = Fz(x);
  return gradient;
}

/** 
 * grad returns the gradient of the function (proc()) for a point or region.
 * Basically, the gradient returns a vector of partial derivatives, or more
 * formally grad(proc()) = (dproc/dx,dproc/dy,dproc/dz).
 * @note Unless grad is overwritten, it defaults to the numerical
 * approximation of the gradient = (proc(x+eps) - proc(x))/eps ...
 * @param x Where to evaluate the gradient of the function.
 * @return The gradient result.
 * @sa setEpsilon(double)
 * @todo remove this function in favor of a template on grad()
 */
Box3d Implicit::grad(const Box<double>& x) 
{ 
  Box3d gradient;
  gradient[0] = Fx(x);
  gradient[1] = Fy(x);
  gradient[2] = Fz(x);
  return gradient;
}

/** 
 * hess returns the Hessian of the function (proc()) for a point or region.
 * Basically, the Hessian is a matrix of second partial derivatives, or more
 * formally hess(proc()) =
 *
 * d^2proc/dx^2, d^2proc/dxdy d^2proc/dxdz \\
 * d^2proc/dxdy, d^2proc/dy^2 d^2proc/dydz \\
 * d^2proc/dxdz, d^2proc/dydz d^2proc/dz^2
 *
 * @note Unless hess is overwritten, it defaults to the numerical
 * approximation of the Hessian.
 * @param x Where to evaluate the Hessian of the function.
 * @return The Hessian result.
 * @sa setEpsilon(double)
 * @note The Hessian is symmetric for our uses.
 */
gmMatrix3 Implicit::hess(const gmVector3 & x) 
{
  gmMatrix3 hessian(Fxx(x),  Fxy(x),  Fxz(x),
                    Fyx(x),  Fyy(x),  Fyz(x),
                    Fzx(x),  Fzy(x),  Fzz(x));
  return hessian;
}

/** 
 * hess returns the Hessian of the function (proc()) for a point or region.
 * Basically, the Hessian is a matrix of second partial derivatives, or more
 * formally hess(proc()) =
 *
 * d^2proc/dx^2, d^2proc/dxdy d^2proc/dxdz \\
 * d^2proc/dxdy, d^2proc/dy^2 d^2proc/dydz \\
 * d^2proc/dxdz, d^2proc/dydz d^2proc/dz^2
 *
 * @note Unless hess is overwritten, it defaults to the numerical
 * approximation of the Hessian.
 * @param x Where to evaluate the Hessian of the function.
 * @return The Hessian result.
 * @sa setEpsilon(double)
 * @note The Hessian is symmetric for our uses.
 * @todo remove this function in favor of a template on hess().
 */
IMatrix3d Implicit::hess(const Box<double>& x) 
{ 
  IMatrix3d hessian(Fxx(x),  Fxy(x),  Fxz(x),
                    Fyx(x),  Fyy(x),  Fyz(x),
                    Fzx(x),  Fzy(x),  Fzz(x));
  return hessian;
}

/**
 * Lipschitz bound of proc. This is an upper bound on the gradient
 * magnitude.  
 * @param x Box over which the Lipschitz bound should be computed. 
 * @note Functions that are not globally Lipschitz may be Lipschitz
 *       over a compact domain. Also a tighter Lipschitz bound might be
 *       used depending on the domain.
 */
double Implicit::lipschitz(const Box<double>& x) 
{
  Box3d v = grad(x);
  double l = gmVector3(v[0].low(),v[1].low(),v[2].low()).lengthSquared();
  double h = gmVector3(v[0].high(),v[1].high(),v[2].high()).lengthSquared();
  return (l > h) ? sqrt(l) : sqrt(h);
}

/** 
 * Deriv of func wrt parameters.
 * Differentiates proc with respect to the parameter vector.
 * @param x Point in space to differentiate proc() at.
 * @param dfdq the returned gradient of proc wrt parameters
 *             (dproc/dq[0],dproc/dq[1],...)
 * @note Default procq is implemented using central differences
 * @sa setEpsilon
 */
void Implicit::procq(const gmVector3 & x, double* dfdq)
{
  // HACK - procq should not upset the parameter "dirty bit" for compact RBFs.
  CompactRBF* crbf = dynamic_cast<CompactRBF*>(this);
	bool pchanged;

  if(crbf)
	 pchanged = crbf->param_changed;

  int len = qlen();
  double* q = new double[len];
  getq(q);
  for (int i = 0; i < len; i++)
    {
      q[i] += m_epsilon; 
      _setq(q);
      dfdq[i] = proc(x);
      q[i] -= 2.0*m_epsilon; 
      _setq(q);
      dfdq[i] -= proc(x);
      dfdq[i] = 1.0 / (2.0*m_epsilon) * dfdq[i];
      q[i] += m_epsilon; 
      _setq(q);
    }
  delete [] q;

	// HACK CONTINUED!
	if(crbf)
		crbf->param_changed = pchanged;
}

///< Evaluation of dFdQ into array
void Implicit::procq(const gmVector3 & x, DoubleArray &result)
{
  procq(x,&result[0]);
}

///< Evaluation of dFdQ into vector
void Implicit::procq(const gmVector3 & x, DoubleVector &result) 
{
  procq(x,&result[0]);
} 

///< Eval of dFdQ into TNT
void Implicit::procq(const gmVector3 & x, TNT::Vector<double> &q) 
{
  std::vector<double> sq(qlen());
  procq(x,sq);
  for (unsigned int i = 0; i < qlen(); i++)
    q[i] = sq[i];
}

void Implicit::procq(const Box<double>& x, Intervald *dfdq)
{
/*
// NOT QUITE SURE HOW TO IMPLEMENT THIS YET.  3D/4D ???
// PROBABLY NEED TO OVERRIDE IN SUBCLASSES TO WORK.
  double *q = (double *) malloc(sizeof(double)*qlen());
  getq(q);
  for (int i = 0; i < qlen(); i++)
    {
      q[i] += m_epsilon; 
      _setq(q);
      dfdq[i] = proc(x);
      q[i] -= 2.0*m_epsilon; 
      _setq(q);
      dfdq[i] -= proc(x);
      dfdq[i] = 1.0 / (2.0*m_epsilon) * dfdq[i];
      q[i] += m_epsilon; 
      _setq(q);
    }
  delete q;
*/
}

void Implicit::gradq(const Box<double>& x, 
                     Intervald *dfdxdq, Intervald *dfdydq, Intervald *dfdzdq)
{
// NO IDEA HOW TO DO THIS IN GENERAL.
// PROBABLY ALWAYS NEED TO OVERRIDE IN SUBCLASSES TO WORK.
}

/**
 * Derivative of the function (proc) wrt t.  In this case, proc is a 4D
 * function (x,y,z,t).  So by using derivation by parts, we get:
 * <pre>
 * dproc(X,t)   dproc(X,t)   dq(t)
 * ---------- = ---------- . -----
 *     dt          dq(t)       dt
 * </pre>
 * dq(t)/dt is dqtdt() which is simply an array of (q_curr[i] - q_old[i]).
 * dproc(X,t)/dq(t) must be defined by each subclass since it is unique for
 * each implicit function/operator.  It would return an array of what is
 * left after you remove the Qs from each element (sort of).  Once you have
 * the two arrays, you take the dot product of them to return a single
 * (Interval) value.
 */
Intervald Implicit::proct(const Box<double>&  b)
{
  Intervald retval(0.0);

  int len = qlen();
  if (len > 0)
    {
      double* dq = new double[len];
      Intervald *pq = (Intervald *) malloc(sizeof(Intervald)*len);
      dqtdt(dq);
      procq(b,pq);

      for (int i = 0; i < len; i++)
        retval += Intervald(dq[i]) * pq[i];

      delete[] dq;
      delete pq;
    }
  return retval;
}

/**
 * Derivative of the gradient of the function (grad) wrt t.  In this case,
 * grad is a 4D function (x,y,z,t).  So by using derivation by parts, we get:
 * <pre>
 * dgrad(X,t)   dgrad(X,t)   dq(t)
 * ---------- = ---------- . -----
 *     dt          dq(t)       dt
 * </pre>
 * dq(t)/dt is dqtdt() which is simply an array of (q_curr[i] - q_old[i]).
 * dgrad(X,t)/dq(t) must be defined by each subclass since it is unique for
 * each implicit function/operator.  It would return an array of what is
 * left after you remove the Qs from each element (sort of).  Once you have
 * the two arrays, you take the dot product of them to return a single
 * (Interval) value.
 */
Box3d Implicit::gradt(const Box<double>&  b)
{
  Box3d retval(0.0);

  int len = qlen();
  if (len > 0)
    {
      Intervald tdq; // Temporary variable
      double* dq = new double[len];
      Intervald *pxq = (Intervald *) malloc(sizeof(Intervald)*len);
      Intervald *pyq = (Intervald *) malloc(sizeof(Intervald)*len);
      Intervald *pzq = (Intervald *) malloc(sizeof(Intervald)*len);
      dqtdt(dq);
      gradq(b,pxq,pyq,pzq);

      for (int i = 0; i < len; i++)
        {
          tdq = Intervald(dq[i]);
          retval[0] += tdq * pxq[i];
          retval[1] += tdq * pyq[i];
          retval[2] += tdq * pzq[i];
        }

      delete[] dq;
      delete pxq;
      delete pyq;
      delete pzq;
    }
  return retval;
}

/** 
 * Derivative of q(t) wrt t.  This is used for 4D critical point finding.
 * The implicit's Qs change from qold[i] (at t = 0) to q[i] (at t = 1).  So,
 * using a LRP, we get q(t)[i] = t * (q[i] - qold[i]) + qold[i].  Taking the
 * derivative wrt t, we are left with q[i] - qold[i].  This is what gets
 * returned in the passed-in array (pointer).  Note that you must allocate
 * memory for dq before calling this method.
 * @param dq The returned array of (q_curr[i] - q_prev[i]).
 * @note This method is called by proct() and gradt().  It's probably not
 * needed by anything else.
 */
void Implicit::dqtdt(double *dq)
{
  int len = qlen();
  if (len > 0)
    {
      double* prevq = new double[len];
      double* currq = new double[len];

      getqold(prevq);
      getq(currq);
      for (int i = 0; i < len; i++)
        dq[i] = currq[i] - prevq[i];

      delete[] prevq;
      delete[] currq;
    }
}

/** 
 * Return the normal of the Implicit function at a point (ie. the 
 * normalized grad()).
 * @param x a point in space
 * @return unit normal at x
 * @note If gradient is zero length, returns (1,0,0).
 */
gmVector3 Implicit::normal(const gmVector3 & x)
{
  gmVector3 n = grad(x);
  
  if (gmIsZero(n.length()))
    n = gmVector3(1.0,0.0,0.0);
  else
    n.normalize();
  
  return n;
}

/** 
 * Returns the tangent of the Implicit function at a point.
 * @param x a point in space
 * @return (0,0,0)
 * @todo implement
 * @todo How is this tangent vector fixed to a consistent direction?
 */
gmVector3 Implicit::tangent(const gmVector3 & x)
{
  gmVector3 tangent(0,0,0);
  return tangent;
}

/** 
 * Returns the binormal of the Implicit function at a point.
 * @param x a point in space
 * @return (0,0,0)
 * @todo implement
 */
gmVector3 Implicit::binormal(const gmVector3 & x)
{
  gmVector3 binormal(0,0,0);
  return binormal;
}

/** 
 * Returns the curvature of the Implicit function at a point.
 * @param x a point in space
 * @return 0
 * @todo implement
 * @todo Which curvature does k return?
 */
double Implicit::k(const gmVector3 & x)
{
  double curv=0;
  return curv;
}

static int usecurvature = 2; // can be reset via the debugger
static bool denominatorfourth = true; // can be reset via the debugger

/** 
 * Returns the Gaussian curvature of the Implicit function at a point.
 * @param x a point in space
 */
double Implicit::gaussianCurvature(const gmVector3 & x)
{
	double fx = Fx(x), fy = Fy(x), fz = Fz(x);

	double den = fx*fx + fy*fy + fz*fz;
	if (denominatorfourth) den *= den;
	if (den == 0.0) return HUGE;	// HUGE defined in math.h
	den = 1.0/den;

	double fxx = Fxx(x), fxy = Fxy(x), fxz = Fxz(x),
						 fyy = Fyy(x), fyz = Fyz(x),
									   fzz = Fzz(x);

	/* this is the version from sethian's level set book, but should probably
	 * have |grad f|^4 in the denominator.
	 */
/*
	double sethian_kg = ((fyy*fzz - fyz*fyz)*fx*fx + (fzz*fxx - fxz*fxz)*fy*fy + (fxx*fyy - fxy*fxy)*fz*fz +
		2.0*fx*fy*(fxz*fyz - fxy*fzz) - 2.0*fy*fz*(fxy*fxz - fyz*fxx) - 2.0*fx*fz*(fxy*fyz - fxz*fyy))
		* den;
*/
	/* This is Spivak's version, taken from Suffern & Balsys 2003
		A = 2[fx fy(fxy fzz - fxz fyz) + fx fz(fxz fyy - fxy fyz) + fy fz(fyz fxx - fxy fxz)]
			- fx^2(fyy fzz - fyz^2) - fy^2(fxx fzz - fxz^2) - fz^2(fxxfyy - fxy^2)
		K = A/(fx^2 - fy^2 + fz^2)^2

		which should be -|hess(x) grad(x)|
		                 |grad(x)   0    |/|grad(x)|^4
	 */

	double spivak_kg = (2.0*(fx*fy*(fxy*fzz - fxz*fyz) + fx*fz*(fxz*fyy - fxy*fyz) + fy*fz*(fyz*fxx - fxy*fxz)) -
					    fx*fx*(fyy*fzz - fyz*fyz) - fy*fy*(fxx*fzz - fxz*fxz) - fz*fz*(fxx*fyy - fxy*fxy))
					   * -den;
					   //-(gradmag2 * gradmag2);

	/* This version is due to Ron Goldman's notes, which update Spivak's derivation. But it doesn't seem
	 * to give the right result.
	 */
/*
	gmMatrix3 hstar = this->hess(x).adjoint();
	gmVector3 grad = this->grad(x);
	//double gradmag4 = grad.lengthSquared();
	//gradmag4 *= gradmag4;

	double ron_kg = dot(grad,hstar * grad) * den;
*/
	//return usecurvature == 1 ? sethian_kg : usecurvature == 2 ? spivak_kg : ron_kg;
	return spivak_kg;
}

Intervald Implicit::gaussianCurvature(const Box<double>& x)
{
	Intervald fx = Fx(x), fy = Fy(x), fz = Fz(x);
	
	Intervald den = fx.squared() + fy.squared() + fz.squared();
	den = den.squared();
	
	Intervald fxx = Fxx(x), fxy = Fxy(x), fxz = Fxz(x),
		fyy = Fyy(x), fyz = Fyz(x),
		fzz = Fzz(x);
	
	return  (2.0*(fx*fy*(fxy*fzz - fxz*fyz) + fx*fz*(fxz*fyy - fxy*fyz) + fy*fz*(fyz*fxx - fxy*fxz)) -
			 fx.squared()*(fyy*fzz - fyz.squared()) - fy.squared()*(fxx*fzz - fxz.squared()) - fz.squared()*(fxx*fyy - fxy.squared())) * (-1.0) / den;
}



double Implicit::numeratorGaussianCurvature(const gmVector3 & x)
{
	double fx = Fx(x), fy = Fy(x), fz = Fz(x);

	double fxx = Fxx(x), fxy = Fxy(x), fxz = Fxz(x),
		fyy = Fyy(x), fyz = Fyz(x),
		fzz = Fzz(x);
	
	return  (2.0*(fx*fy*(fxy*fzz - fxz*fyz) + fx*fz*(fxz*fyy - fxy*fyz) + fy*fz*(fyz*fxx - fxy*fxz)) -
					    fx*fx*(fyy*fzz - fyz*fyz) - fy*fy*(fxx*fzz - fxz*fxz) - fz*fz*(fxx*fyy - fxy*fxy));
}

Intervald Implicit::numeratorGaussianCurvature(const Box<double>& x)
{
	Intervald fx = Fx(x), fy = Fy(x), fz = Fz(x);
	
	Intervald fxx = Fxx(x), fxy = Fxy(x), fxz = Fxz(x),
		fyy = Fyy(x), fyz = Fyz(x),
		fzz = Fzz(x);
	
	return   (2.0*(fx*fy*(fxy*fzz - fxz*fyz) + fx*fz*(fxz*fyy - fxy*fyz) + fy*fz*(fyz*fxx - fxy*fxz)) -
			  fx.squared()*(fyy*fzz - fyz.squared()) - fy.squared()*(fxx*fzz - fxz.squared()) - fz.squared()*(fxx*fyy - fxy.squared()));
}


/*! This method approximates the gradient of gaussian curvature.
 *  \param x A point in space.
 *  \return An approximation of the gradient of gaussian curvature at x.
 *  \author Jared Hoberock
 */
gmVector3 Implicit::gradGaussianCurvature(const gmVector3 & x)
{
  gmVector3 result(0,0,0);

  double step = getEpsilon();

  result[0] = gaussianCurvature(x + gmVector3(step,0,0)) - gaussianCurvature(x);
  result[1] = gaussianCurvature(x + gmVector3(0,step,0)) - gaussianCurvature(x);
  result[2] = gaussianCurvature(x + gmVector3(0,0,step)) - gaussianCurvature(x);
  result /= step;

  return result;
} // end Implicit::gradGaussianCurvature()

/** 
 * Returns the mean curvature of the Implicit function at a point.
 * @param x a point in space
 */
double Implicit::meanCurvature(const gmVector3 & x)
{
	double fx = Fx(x), fy = Fy(x), fz = Fz(x);

	double fxx = Fxx(x), fxy = Fxy(x), fxz = Fxz(x),
						 fyy = Fyy(x), fyz = Fyz(x),
									   fzz = Fzz(x);

	double km = 0.5*((fyy+fzz)*fx*fx + (fzz+fxx)*fy*fy + (fxx+fyy)*fz*fz -
		2.0*fx*fy*fxy - 2.0*fy*fz*fyz - 2.0*fz*fx*fxz)/pow(fx*fx + fy*fy + fz*fz,1.5);

	return km;
}

/*! This method approximates the gradient of mean curvature.
 *  \param x A point in space.
 *  \return An approximation of the gradient of mean curvature at x.
 *  \author Jared Hoberock
 */
gmVector3 Implicit::gradMeanCurvature(const gmVector3 & x)
{
  gmVector3 result(0,0,0);

  double step = getEpsilon();

  result[0] = meanCurvature(x + gmVector3(step,0,0)) - meanCurvature(x);
  result[1] = meanCurvature(x + gmVector3(0,step,0)) - meanCurvature(x);
  result[2] = meanCurvature(x + gmVector3(0,0,step)) - meanCurvature(x);
  result /= step;

  return result;
} // end Implicit::gradMeanCurvature()

/*! This method applies the Shape Operator of the Implicit at point x in direction v.
 *  \param x A point in space.
 *  \param v A vector (not necessarily unit length) to evaluate with respect to.
 *  \see http://mathworld.wolfram.com/ShapeOperator.html
 *  \return A gmVector3, the Shape Operator at x wrt v.
 *
 *  The Shape Operator of a surface M in direction v is defined as
 *  $S(v) = -Dv N$, where $N$ is the unit normal vector field of M, and Dv differentiates in
 *  direction v.  To apply the Shape Operator at point x, in direction v, we evaluate
 *  Dv N by taking the Hessian at x and dotting it with v to yield a three vector.
 *  \author Jared Hoberock
 */
gmVector3 Implicit::shapeOperator(const gmVector3 & x, const gmVector3 & v)
{
  // get the Hessian at x
  gmMatrix3 hessian = hess(x);

  // return the negation of hessian . v
  return - hessian * v;
} // end Implicit::shapeOperator()

/*! 
 *  This method returns the normal curvature at x in direction w.
 *  \param x A point in space.
 *  \param v A vector tangent to the surface (not necessarily of unit length) to evaluate with respect to.
 *  \see http://mathworld.wolfram.com/NormalCurvature.html
 *  \return A scalar, the normal curvature at x in direction w.
 *
 *  Mathworld tells us normal curvature wrt v is the Shape Operator S(v) applied to
 *  the negative derivative of N in direction w (-Dv N).  
 */
double Implicit::normalCurvature(const gmVector3 & x, const gmVector3 & v)
{
  return dot(shapeOperator(x,v), v);
} // end Implicit::radialCurvature()

double Implicit::radialCurvature(const gmVector3 & x, const gmVector3 & v)
{
	gmVector3 normal = this->normal(x);
	gmVector3 w = v - dot(normal, v) * normal;
	double l2w = w.lengthSquared();
	if(l2w == 0) return HUGE;
	return  (1.0/l2w) * dot(w,this->hess(x) * w);
}

double Implicit::numeratorOfRadialCurvature(const gmVector3 & x, const gmVector3 & v)
{
	gmVector3 normal = this->normal(x);
	gmVector3 w = v - dot(normal, v) * normal;
	return  dot(w,this->hess(x) * w);
}

void Implicit::gradRadialCurvature(gmVector3& gradRadCurv, const gmVector3 & x, const gmVector3 & v)
{
	gmMatrix3 Hess = this->hess(x) ; //hessian
	gmVector3 hessv = Hess * v;
	gmVector3 normal = this->normal(x);
	//the projection of the view vector onto the tangent plane
	gmVector3 w = v - dot(normal, v) * normal;
	gmVector3 hessw = Hess * w;
	gmVector3 gradient = this->grad(x);
	double l2w = w.lengthSquared();
	double l2grad = gradient.lengthSquared();
	double l2winv = 1.0/l2w;
	double l2gradinv = 1.0/l2grad;
	double cf =	dot(w,hessw) * 2.0 * l2winv * l2winv;
	double gradfw = dot(gradient, w);
	double gradfv = dot(gradient, v);
	gmVector3 firstt = l2gradinv * (gradfw * (hessv - gradient) + gradfv * hessw);
	gmVector3 secondt = -2 * l2gradinv * l2gradinv *  gradfv * gradfw * (Hess * gradient);
	gmVector3 first = cf * (w + firstt + secondt);
	gmVector3 second(
						 dot((this->hessi(x,0) * w),w),
						 dot((this->hessi(x,1) * w),w),
						 dot((this->hessi(x,2) * w),w)
						 );
		
	gmVector3 third = - hessw - l2gradinv * (dot(hessw, gradient) * (hessv - gradient) 
		+ gradfv * Hess * hessw) + 
		2 * l2gradinv * l2gradinv * (gradfv * dot(hessw,gradient) * hessv);
	
	//see the formula in the briliant paper: "Clip Art Rendering of Smooth Isosurfaces"
	gradRadCurv = first + l2winv * (second + 2 * third); 
		

}

void Implicit::gradNumeratorOfRadialCurvature(gmVector3& gradRadCurv, const gmVector3 & x, const gmVector3 & v)
{
	gmMatrix3 Hess = this->hess(x) ; //hessian
	gmVector3 hessv = Hess * v;
	gmVector3 normal = this->normal(x);
	//the projection of the view vector onto the tangent plane
	gmVector3 w = v - dot(normal, v) * normal;
	gmVector3 hessw = Hess * w;
	gmVector3 gradient = this->grad(x);
	double l2grad = gradient.lengthSquared();
	double l2gradinv = 1.0/l2grad;
	double gradfv = dot(gradient, v);
	gmVector3 second(
						 dot((this->hessi(x,0) * w),w),
						 dot((this->hessi(x,1) * w),w),
						 dot((this->hessi(x,2) * w),w)
						 );
		
	gmVector3 third = - hessw - l2gradinv * (dot(hessw, gradient) * (hessv - gradient) 
		+ gradfv * Hess * hessw) + 
		2 * l2gradinv * l2gradinv * (gradfv * dot(hessw,gradient) * hessv);
	
	//see the formula in the briliant paper: "Clip Art Rendering of Smooth Isosurfaces"
	gradRadCurv = second + 2 * third; 
		

}


/*! This method approximates the gradient of normal curvature.
 *  \param x A point in space.
 *  \param v A vector in the tangent space of this Implicit.
 *  \return An approximation of the gradient of normal curvature at x in direction v.
 *  \author Jared Hoberock
 *
 *  First, normal curvature is defined as: \n
 *  kr(x,v) = (-Hess f(x) . v) . v  \n
 *  We want to evaluate the gradient of this: \n
 *  grad kr(x,v) = grad ((-Hess f(x) . v) . v) \n
 *  grad kr(x,v) = [-grad Hess f(x) . v - Hess f(x) . grad v] . v + [-Hess f(x) . v] . grad v \n
 *  \n
 *  v is a vector in the tangent space of the surface "pointing towards" the camera.  If we assume the camera
 *  is far away from the surface, then v shouldn't change much, and grad v -> 0.  I'm not sure if this is really
 *  correct, since v doesn't actually point towards the camera, it is the view vector projected into the tangent
 *  space of S, but whatever.  So what we have now is:
 *  \n
 *  grad kr(x,v) = -grad Hess f(x) . v . v \n
 *  \n
 *  So what we need to do is approximate the gradient of the Hessian at x and then do a couple of dot products.
 */
gmVector3 Implicit::gradNormalCurvature(const gmVector3 & x, const gmVector3 & v)
{
  // The gradient of the Hessian is a 3x3x3 matrix.  We don't really want to construct this, so we'll construct the matrix
  // grad Hess f(x) . v instead.

  // grad Hess . v = | Hess_x . v |
  //                 | Hess_y . v |
  //                 | Hess_z . v |

  gmMatrix3 hess_x = hessi(x,0);
  gmMatrix3 hess_y = hessi(x,1);
  gmMatrix3 hess_z = hessi(x,2);

  // now dot v against these three matrices to create grad hess . v
  gmVector3 hess_xDotV = hess_x * v;
  gmVector3 hess_yDotV = hess_y * v;
  gmVector3 hess_zDotV = hess_z * v;

  gmMatrix3 gradHessDotV(hess_xDotV[0], hess_xDotV[1], hess_xDotV[2],
                         hess_yDotV[0], hess_yDotV[1], hess_yDotV[2],
                         hess_zDotV[0], hess_zDotV[1], hess_zDotV[2]);

  // the final result is - gradHessDotV . v
  return - gradHessDotV * v;
} // end Implicit::gradNormalCurvature()

/** 
 * Sets the epsilon used by the numerical differentiation methods for grad()
 * and hess().
 * @param epsilon Size of neighborhood for numerical differentiation
 */
void Implicit::setEpsilon(double epsilon)
{
  m_epsilon = epsilon;
}

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %  First Partial Derivative functions                        %    
 %  TODO later replace with better approximations             %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
double Implicit::Fx(const gmVector3 & x) 
{
  gmVector3 temp(x);
  temp[0] += m_epsilon;
  return (proc(temp)-proc(x))/m_epsilon;
}
double Implicit::Fy(const gmVector3 & x)
{
  gmVector3 temp(x);
  temp[1] += m_epsilon;
  return (proc(temp)-proc(x))/m_epsilon;
}
double Implicit::Fz(const gmVector3 & x)
{
  gmVector3 temp(x);
  temp[2] += m_epsilon;
  return (proc(temp)-proc(x))/m_epsilon;
}
Intervald Implicit::Fx(const Box<double>& x) 
{
  Box3d temp(x);
  temp[0] += m_epsilon;
  return (proc(temp)-proc(x))/Intervald(m_epsilon);
}
Intervald Implicit::Fy(const Box<double>& x)
{
  Box3d temp(x);
  temp[1] += m_epsilon;
  return (proc(temp)-proc(x))/Intervald(m_epsilon);
}
Intervald Implicit::Fz(const Box<double>& x)
{
  Box3d temp(x);
  temp[2] += m_epsilon;
  return (proc(temp)-proc(x))/Intervald(m_epsilon);
}

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %  Second Partial Derivative functions wrt X                 %    
 %  TODO later replace with better approximations             %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
double Implicit::Fxx(const gmVector3 & x)
{
  gmVector3 tempXplusE(x),tempXminusE(x);
  tempXplusE[0] += m_epsilon;
  tempXminusE[0] -= m_epsilon;
  return ((proc(tempXminusE)-2*proc(x)+proc(tempXplusE))/(m_epsilon*m_epsilon));
}
double Implicit::Fxy(const gmVector3 & x)
{
  gmVector3 temp(x);
  temp[1] += m_epsilon;
  return (Fx(temp)-Fx(x))/m_epsilon;
}
double Implicit::Fxz(const gmVector3 & x)
{
  gmVector3 temp(x);
  temp[2] += m_epsilon;
  return (Fx(temp)-Fx(x))/m_epsilon;
}
Intervald Implicit::Fxx(const Box<double>& x)
{
  Box3d tempXplusE(x);
  Box3d tempXminusE(x);
  tempXplusE[0] += m_epsilon;
  tempXminusE[0] -= m_epsilon;
  return ((proc(tempXminusE)-Intervald(2.0)*proc(x)+
           proc(tempXplusE))/Intervald(m_epsilon*m_epsilon));
}
Intervald Implicit::Fxy(const Box<double>& x)
{
  Box3d temp(x);
  temp[1] += m_epsilon;
  return (Fx(temp)-Fx(x))/Intervald(m_epsilon);
}
Intervald Implicit::Fxz(const Box<double>& x)
{
  Box3d temp(x);
  temp[2] += m_epsilon;
  return (Fx(temp)-Fx(x))/Intervald(m_epsilon);
}

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %  Second Partial Derivative functions wrt Y                 %    
 %  TODO later replace with better approximations             %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
double Implicit::Fyx(const gmVector3 & x)
{
  gmVector3 temp(x);
  temp[0] += m_epsilon;
  return (Fy(temp)-Fy(x))/m_epsilon;
}
double Implicit::Fyy(const gmVector3 & x)
{
  gmVector3 tempXplusE(x),tempXminusE(x);
  tempXplusE[1] += m_epsilon;
  tempXminusE[1] -= m_epsilon;
  return ((proc(tempXminusE)-2*proc(x)+proc(tempXplusE))/(m_epsilon*m_epsilon));
}
double Implicit::Fyz(const gmVector3 & x)
{
  gmVector3 temp(x);
  temp[2] += m_epsilon;
  return (Fy(temp)-Fy(x))/m_epsilon;
}
Intervald Implicit::Fyx(const Box<double>& x)
{
  Box3d temp(x);
  temp[0] += m_epsilon;
  return (Fy(temp)-Fy(x))/Intervald(m_epsilon);
}
Intervald Implicit::Fyy(const Box<double>& x)
{
  Box3d tempXplusE(x);
  Box3d tempXminusE(x);
  tempXplusE[1] += m_epsilon;
  tempXminusE[1] -= m_epsilon;
  return ((proc(tempXminusE)-Intervald(2.0)*proc(x)+
           proc(tempXplusE))/Intervald(m_epsilon*m_epsilon));
}
Intervald Implicit::Fyz(const Box<double>& x)
{
  Box3d temp(x);
  temp[2] += m_epsilon;
  return (Fy(temp)-Fy(x))/Intervald(m_epsilon);
}

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %  Second Partial Derivative functions wrt Z                 %    
 %  TODO later replace with better approximations             %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
double Implicit::Fzx(const gmVector3 & x)
{
  gmVector3 temp(x);
  temp[0] += m_epsilon;
  return (Fz(temp)-Fz(x))/m_epsilon;
}
double Implicit::Fzy(const gmVector3 & x)
{
  gmVector3 temp(x);
  temp[1] += m_epsilon;
  return (Fz(temp)-Fz(x))/m_epsilon;
}
double Implicit::Fzz(const gmVector3 & x)
{
  gmVector3 tempXplusE(x),tempXminusE(x);
  tempXplusE[2] += m_epsilon;
  tempXminusE[2] -= m_epsilon;
  return ((proc(tempXminusE)-2*proc(x)+proc(tempXplusE))/(m_epsilon*m_epsilon));
}
Intervald Implicit::Fzx(const Box<double>& x)
{
  Box3d temp(x);
  temp[0] += m_epsilon;
  return (Fz(temp)-Fz(x))/Intervald(m_epsilon);
}
Intervald Implicit::Fzy(const Box<double>& x)
{
  Box3d temp(x);
  temp[1] += m_epsilon;
  return (Fz(temp)-Fz(x))/Intervald(m_epsilon);
}
Intervald Implicit::Fzz(const Box<double>& x)
{
  Box3d tempXplusE(x);
  Box3d tempXminusE(x);
  tempXplusE[2] += m_epsilon;
  tempXminusE[2] -= m_epsilon;
  return ((proc(tempXminusE)-Intervald(2.0)*proc(x)+
           proc(tempXplusE))/Intervald(m_epsilon*m_epsilon));
}

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %  Third Partial Derivative functions wrt X                 %    
 %  TODO later replace with better approximations             %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
double Implicit::Fxxi(const gmVector3 & x, unsigned int i)
{
  gmVector3 temp(x);
  temp[i] += m_epsilon;
  return (Fxx(temp)-Fxx(x))/m_epsilon;
} // end Implicit::Fxxi()

double Implicit::Fxyi(const gmVector3 & x, unsigned int i)
{
  gmVector3 temp(x);
  temp[i] += m_epsilon;
  return (Fxy(temp) - Fxy(x)) / m_epsilon;
} // end Implicit::Fxyi()

double Implicit::Fxzi(const gmVector3 & x, unsigned int i)
{
  gmVector3 temp(x);
  temp[i] += m_epsilon;
  return (Fxz(temp) - Fxz(x)) / m_epsilon;
} // end Implicit::Fxzi()

double Implicit::Fyxi(const gmVector3 & x, unsigned int i)
{
  gmVector3 temp(x);
  temp[i] += m_epsilon;
  return (Fyx(temp) - Fyx(x)) / m_epsilon;
} // end Implicit::Fyxi()

double Implicit::Fyyi(const gmVector3 & x, unsigned int i)
{
  gmVector3 temp(x);
  temp[i] += m_epsilon;
  return (Fyy(temp) - Fyy(x)) / m_epsilon;
} // end Implicit::Fyyi()

double Implicit::Fyzi(const gmVector3 & x, unsigned int i)
{
  gmVector3 temp(x);
  temp[i] += m_epsilon;
  return (Fyz(temp) - Fyz(x)) / m_epsilon;
} // end Implicit::Fyzi()

double Implicit::Fzxi(const gmVector3 & x, unsigned int i)
{
  gmVector3 temp(x);
  temp[i] += m_epsilon;
  return (Fzx(temp) - Fzx(x)) / m_epsilon;
} // end Implicit::Fzxi()

double Implicit::Fzyi(const gmVector3 & x, unsigned int i)
{
  gmVector3 temp(x);
  temp[i] += m_epsilon;
  return (Fzy(temp) - Fzy(x)) / m_epsilon;
} // end Implicit::Fzyi()

double Implicit::Fzzi(const gmVector3 & x, unsigned int i)
{
  gmVector3 temp(x);
  temp[i] += m_epsilon;
  return (Fzz(temp) - Fzz(x)) / m_epsilon;
} // end Implicit::Fzzi()

//interval versions
Intervald Implicit::Fxxi(const Box<double>& x, unsigned int i)
{
	Box3d temp(x);
	temp[i] += m_epsilon;
	return (Fxx(temp)-Fxx(x))/Intervald(m_epsilon);
} // end Implicit::Fxxi()

Intervald Implicit::Fxyi(const Box<double>& x, unsigned int i)
{
	Box3d temp(x);
	temp[i] += m_epsilon;
	return (Fxy(temp) - Fxy(x)) / Intervald(m_epsilon);
} // end Implicit::Fxyi()

Intervald Implicit::Fxzi(const Box<double>& x, unsigned int i)
{
	Box3d temp(x);
	temp[i] += m_epsilon;
	return (Fxz(temp) - Fxz(x)) / Intervald(m_epsilon);
} // end Implicit::Fxzi()

Intervald Implicit::Fyxi(const Box<double>& x, unsigned int i)
{
	Box3d temp(x);
	temp[i] += m_epsilon;
	return (Fyx(temp) - Fyx(x)) / Intervald(m_epsilon);
} // end Implicit::Fyxi()

Intervald Implicit::Fyyi(const Box<double>& x, unsigned int i)
{
	Box3d temp(x);
	temp[i] += m_epsilon;
	return (Fyy(temp) - Fyy(x)) / Intervald(m_epsilon);
} // end Implicit::Fyyi()

Intervald Implicit::Fyzi(const Box<double>& x, unsigned int i)
{
	Box3d temp(x);
	temp[i] += m_epsilon;
	return (Fyz(temp) - Fyz(x)) / Intervald(m_epsilon);
} // end Implicit::Fyzi()

Intervald Implicit::Fzxi(const Box<double>& x, unsigned int i)
{
	Box3d temp(x);
	temp[i] += m_epsilon;
	return (Fzx(temp) - Fzx(x)) / Intervald(m_epsilon);
} // end Implicit::Fzxi()

Intervald Implicit::Fzyi(const Box<double>& x, unsigned int i)
{
	Box3d temp(x);
	temp[i] += m_epsilon;
	return (Fzy(temp) - Fzy(x)) / Intervald(m_epsilon);
} // end Implicit::Fzyi()

Intervald Implicit::Fzzi(const Box<double>& x, unsigned int i)
{
	Box3d temp(x);
	temp[i] += m_epsilon;
	return (Fzz(temp) - Fzz(x)) / Intervald(m_epsilon);
} // end Implicit::Fzzi()

////////////////////////
Intervald Implicit::D3F(const Box<double>& x, unsigned int i, unsigned int j, unsigned int k)
{
	Box3d temp(x);
	temp[i] += m_epsilon;
	return (hess(temp)(j+1,k+1) - hess(x)(j+1,k+1)) / Intervald(m_epsilon);
} // end Implicit::D3F()
/////////////////////////////
double Implicit::D3F(const gmVector3& x,  unsigned int i,  unsigned int j,  unsigned int k) 
{
	gmVector3 temp(x);
	temp[i] += m_epsilon;
	return (hess(temp)[j][k] - hess(x)[j][k]) / m_epsilon;
} // end Implicit::D3F()
/////////////////////////////

gmMatrix3 Implicit::hessi(const gmVector3 &x, const unsigned int i)
{
  gmVector3 temp(x);
  temp[i] += m_epsilon;
  return (hess(temp) - hess(x)) / m_epsilon;
} // end Implicit::hessi()

//see comment in h file.
/// the gradient of the Norm Squared (gradient)
//gmVector3 Implicit::gradGradN2(const gmVector3 &x)
//{
//	double fx=Fx(x);
//	double fy=Fy(x);
//	double fz=Fz(x);
//
//	// Fxz == Fzx
//	double fxy=Fxy(x);
//	double fxz=Fxz(x);
//	double fyz=Fyz(x);
//	
//    return 2*gmVector3(
//		fx*Fxx(x)+fy*fxy+fz*fxz,
//		fx*fxy+fy*Fyy(x)+fz*fyz,
//		fx*fxz+fy*fyz+fz*Fzz(x));
//}


//
///// Polygonizes the Implicit into the mesh stored in Surface::m_mesh
//char *Implicit::polygonize(double size, int bounds, gmVector3 x, 
//                           int (*triproc)(int,int,int,jbVertices))
//{
//  jbProcess p;
//  return p.polygonize(this, size, bounds, x, triproc);
//}
//
///// Polygonizes specifically using marching cubes
//char *Implicit::marchingcubes(double size, gmVector3 x0, gmVector3 x1, 
//                              int (*triproc)(int,int,int,jbVertices))
//{
//  jbProcess p;
//  return p.marchingcubes(this, size, x0, x1, triproc);
//}

///# of implicit model parameters for this object
unsigned int Implicit::plen()
{  
  int sublen = 0;
  std::vector<Implicit*> children;
  
  getChildren(children);
  for (unsigned int i = 0; i < children.size(); i++)
    {
//      Implicit* temp = children[i];
      if (children[i])
        sublen += children[i]->qlen();    
    }
  
  return qlen() - sublen;
}

/// # of implicit model parameters for the tree
unsigned int Implicit::qlen()
{
  return 0;
}

/** The name of each parameter.  The names of the parameters of the function
 * (if any) are put into an array which gets returned.
 * @param qn An array of size qlen of strings listing names of parameters
 */
void Implicit::getqname(char**)
{ 
}

/**
 * Fetches parameter names into a vector.
 * @param names Vector of names to update.
 */
void Implicit::getqname(NameVector &names) 
{
  names.clear();
  
  int theQlen = qlen();
  names.reserve(theQlen);
  char** n = new char*[theQlen];
  
  getqname(n);
  
  for(int i=0;i<theQlen;i++) 
    names.push_back(n[i]);

  delete [] n;
}


/**
 * Solves a system of equations using Cholesky factorization.
 * @param A Matrix of coefficients.
 * @param x Vector of variables.
 * @param b Vector of values.
 * @return True if successful.
 */
bool solveCholesky(TNT::Matrix<double>& A, TNT::Vector<double>& x, TNT::Vector<double>& b) 
{
	int n = x.size();
	TNT::Matrix<double> L(n, n);
	if (TNT::Cholesky_upper_factorization(A, L) != 0) return false;
	TNT::Vector<double> y = TNT::Lower_triangular_solve(L, b);
	TNT::Transpose_View<TNT::Matrix<double> > T = TNT::Transpose_view(L);
	x = TNT::Upper_triangular_solve(T, y);
	return true;
}

void Implicit::interpolate(Particles *p, float phi)
{
	unsigned int i,j,k;
	gmVector3 gradient;
	gmVector3 xi;
	unsigned int qlen = this->qlen();	// # of parameters
	unsigned int n_free;					// # of free parameters
	unsigned int n = p->size();			// # of control particles

	ParticlePosition *position;
	ParticleVelocity *velocity;

	// This assumes Particles p has only one position attribute
	p->findAttribute(position);
	if (!position) return;

	p->findAttribute(velocity);
	if (!velocity) return;

	// for now just make everything flexible
	std::valarray<bool> flexible(true,qlen);

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
		procq(position->getPosition(i), dFdQi);

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
		//b[i] += dot(grad(position->getPosition(i)), velocity->v[i]) + 15 * proc(position->getPosition(i));
		b[i] += dot(grad(position->getPosition(i)), velocity->v[i]) + phi * proc(position->getPosition(i));

		// Build row i of matrix A
		for(j=0; j<n; j++)
		{
			A(i+1,j+1) = 0.0;

			// Get derivative of F wrt Q at particle j
			procq(position->getPosition(j), dFdQj);

			// Only the subset of this array corresponding to flexible
			// parameters is used.
			for (k = 0; k < qlen; k++)
			if (flexible[k]) 
			A(i+1,j+1) += dFdQi[k]*dFdQj[k];
		} // end matrix loop
	} // end control particle loop

	#ifdef DEBUG_MATRIX
	std::cerr << A;
	#endif

	// Try solving Ax = b using Cholesky
	if (!solveCholesky(A, x, b))
	{
		SVD svd;
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
		procq(position->getPosition(j), dFdQj);
		dFdQj *= (double)x[j];

		// Only modify those that are flexible
		for(k=0; k<qlen; k++) 
			if (flexible[k])
				dqdt[k] -= dFdQj[k];
	}

	// Apply the changes to the implicit surface parameter
	std::valarray<double> q(qlen);
	getq(q);
	q += dqdt * (double)p->dt;
	setq(q);
}


/// Number of children currently defined.
int Implicit::numChildren() 
{ 
  return 0; 
}

/** Max number of children that can be defined.  This will return -1 for
 * operators that can have an arbitrary number of operands.
 */
int Implicit::maxChildren() 
{ 
  return 0; 
}

/**
 * Get copy of the implicit model parameters.
 * @param A pointer to the array to be returned.
 */
void Implicit::getq(double* q)
{
}

/**
 * Get parameters into a double array.
 * @param q DoubleArray for parameters.
 */
void Implicit::getq(DoubleArray &q) 
{
  if (q.size()!=qlen()) 
    q.resize(qlen());
  
  getq(&q[0]);
} 

/**
 * Get parameters into a double vector.
 * @param q DoubleVector for parameters.
 */
void Implicit::getq(DoubleVector &q) 
{
  unsigned int theQlen = qlen();
  q.reserve(theQlen);
  if (q.size()!=theQlen) 
    q.resize(theQlen);
  
  getq(&q[0]);
}

/**
 * Get parameters into TNT double vector.
 * @param q TNT::Vector<double> for parameters.
 */
void Implicit::getq(TNT::Vector<double> &q)
{
  std::vector<double> sq(qlen());
  getq(sq);
  for (unsigned int i = 0; i < qlen(); i++)
    q[i] = sq[i];
}

/**
 * Get the OLD implicit model parameters.
 * @param A pointer to the array to be returned.
 */
void Implicit::getqold(double* q)
{
  for (unsigned int i = 0; i < qlen(); i++)
    q[i] = qold[i];
}

/**
 * Get OLD parameters into a double array.
 * @param q DoubleArray for parameters.
 */
void Implicit::getqold(DoubleArray &q) 
{
  if (q.size()!=qlen()) 
    q.resize(qlen());
  
  getqold(&q[0]);
} 

/**
 * Get OLD parameters into a double vector.
 * @param q DoubleVector for parameters.
 */
void Implicit::getqold(DoubleVector &q) 
{
  unsigned int theQlen = qlen();
  q.reserve(theQlen);
  if (q.size()!=theQlen) q.resize(theQlen);
  
  getqold(&q[0]);
}

/**
 * Get OLD parameters into TNT double vector.
 * @param q TNT::Vector<double> for parameters.
 */
void Implicit::getqold(TNT::Vector<double> &q)
{
  std::vector<double> sq(qlen());
  getqold(sq);
  for (unsigned int i = 0; i < qlen(); i++)
    q[i] = sq[i];
}

/** 
 * Is qold different from current Qs?  Basically, we get the current Qs and
 * the previous Qs (qold) and compare them to see if any of the Qs have
 * changed.  If so, we return true.
 * @return True of the current Qs are different from the previous Qs.
 */
#if 0
bool Implicit::qsChanged(void)
{
  bool retval = false;
  int len = qlen();
  if (len > 0)
    {
      double* prevq = new double[len];
      double* currq = new double[len];

      getqold(prevq);
      getq(currq);
      int i = 0;
      while ((i < len) && (retval == false))
        {
          if (prevq[i] != currq[i])
            retval = true;
          i++;
        }
      
      delete[] prevq;
      delete[] currq;
    }
  return retval;
}
#else
/** 
 * Is qold different from current Qs?  This is a slightly stupid function
 * that checks to see if the qsChangedFlag was set to "true" by setq().  If
 * so, then it returns true and resets the qsChangedFlag to "false".  Thus,
 * if you call qsChanged multiple times before you call setq again,
 * qsChanged will return "true" only once.  Note that this incorrectly
 * returns "true" if setq() was called with the same set of Qs that were
 * previously stored.  Hopefully that won't occur too often.
 * @return True of the current Qs are different from the previous Qs.
 */
bool Implicit::qsChanged(void)
{
  bool retval = qsChangedFlag;
  qsChangedFlag = false;
  return retval;
}
#endif

/** Get memory space for old "array".  This method is called by
 *  setq(double*) to allocate memory for qold AND set the value of
 *  qoldsize.  This is necessary because some of the implicits change the
 *  size of their Qs during runtime.  For example, Intersection initially
 *  has 0 Qs.  When one child is set, then it has a qlen of that child.
 *  When both children are set, it has a qlen of the sum of the qlens of the
 *  children.  We need to reallocate memory for qold when this happens.
 */
void Implicit::allocateQold(void)
{
  if (qold == NULL)
    {
      qoldsize = qlen();  // Keep track of current size of qold array
      qold = new double[qoldsize];
    }
 
  // Check to see if we allocated enough space for qold based on the 
  // current qlen().  If not, free up qold and allocate new memory.
  if ((qold != NULL) && (qoldsize != qlen()))
    { 
      delete[] qold;
      qoldsize = qlen();
      qold = new double[qoldsize];
    }
}

/** Change the implicit model parameter vector.  This method first saves the
 *  last Qs into qold and then calls the object's _setq() to set the current
 *  Q array.  If you don't need to save the last Qs, then override this
 *  method.  Otherwise, override _setq() which does the actual work of
 *  setting the Q array and doesn't save qold.
 *  @param A pointer to the array holding the new Qs
 *  @see double* qold
 *  @see _setq(double*)
 */
void Implicit::setq(double* q)
{
  // Make sure enough space is allocated for qold.
  allocateQold();
  // Always set the qold before setting the current Qs.
  if (qold != NULL)
    getq(qold);
  qsChangedFlag = true;
  _setq(q);
}

/**
 * Set parameters from a double array.
 * @param q DoubleArray of parameters.
 */
void Implicit::setq(DoubleArray &q)
{  
  if (q.size() != qlen())
    q.resize(qlen(), 0.0);
  
  setq(&q[0]);
}

/**
 * Set parameters from a double vector.
 * @param q DoubleVector of parameters.
 */
void Implicit::setq(DoubleVector &q)
{  
  if (q.size() != qlen())
    q.resize(qlen(), 0.0);
  
  setq(&q[0]);
}

/**
 * Set parameters from a double vector.
 * @param q TNT::Vector<double> for parameters.
 */
void Implicit::setq(TNT::Vector<double> &q)
{
  DoubleVector sq(qlen());
  for (unsigned int i = 0; i < qlen(); i++)
    sq[i] = q[i];
  setq(sq);
}

/** Overridden by subclasses to set Qs.  This is the main method that needs
 *  to be defined in all of Implicit's subclasses.  It does the actual work
 *  of setting the Q parameters for the Implicit.  This is called by the
 *  main Implicit::setq(double*) which saves the qold array first and then
 *  calls the object's _setq to set the current parameter array.
 *  @param A pointer to the array holding the new Qs
 *  @note Subclasses should override _setq() to set the actual Q array
 *  @see setq(double*)
 */
void Implicit::_setq(double* q)
{
}

/**
 * If this implicit has children, this function sets them.
 * @param  children The children of this object.
 * @return False if the child vector is invalid.
 */
bool Implicit::setChildren(std::vector<Implicit*> children)
{
  if (children.size() != maxChildren())
    return false;

  int max = (maxChildren() == -1) ? numChildren() : maxChildren();

  for (int i = 0; i < max; i++)
    setChild(i, children[i]);
  
  return true;
}

/**
 * If this implicit has children, this function sets the child specified by
 * the index.
 * @param  index The index of the child to set.
 * @param  child The child of this object.
 * @return False if the child is invalid.
 * @note    Needs to be overridden.
 */
bool Implicit::setChild(int index, Implicit* child) 
{ 
  return false; 
}

/**
 * Get the specified child.  
 * @param  index The index of the child to return
 * @return NULL if no children or child[index]
 * @note    Needs to be overridden.
 */
Implicit* Implicit::getChild(int index) 
{ 
  return NULL; 
}

/** 
 * Retrieve the children of this implicit.
 * This is a convenience function that will return all of the Implicit's
 * children into a vector. Subclasses shouldn't need to overwrite this
 * function as it simply uses getChild.
 * @param children  The vector to put all the children into.
 * @sa getChild(int)
 */
void Implicit::getChildren(std::vector<Implicit*> &children)
{
  for (int i = 0; i < maxChildren(); i++)
    children.push_back(getChild(i));
}

/**
 * Returns a char * in the XPM image format. This is a icon that can be used
 * to represent the implicit in a GUI environment.
 */
const char ** Implicit::getPixmapXPM(const int& size) const
{
  if (size <= 16)
    return (const char **)default_pixmap16;
  else if (size <= 32)
    return (const char **)default_pixmap32;
  else
    return (const char **)default_pixmap48;
}

void Implicit::output(std::ostream &out)
{
	/* output only first plen of the qlen q-parameters */
	int qlen = this->qlen();
	int plen = this->plen();
	std::valarray<double> q(qlen);
	getq(q);
	for (int i = 0; i < plen; i++) {
		out << (i ? ' ' : '\t') << q[i];
	}
	if (plen)
		out << "\n";
}

void Implicit::input(std::istream &in)
{
	/* Only read in first plen of the qlen q-parameters */
	int plen = this->plen();
	int qlen = this->qlen();
	DoubleVector q(qlen);
	getq(q);
	for (int i = 0; i < plen; i++)
		in >> q[i];
	setq(q);
}

/** This is defined here in Implicit.cpp because we need to construct an Implicit by name
 * and need the implicit factory for that.
 */
std::istream &operator>>(std::istream &in, Surfaces *surfaces)
{
	std::string surfclass, surfname;
	char c;
	ImplicitRegistry *impReg = &IMPLICIT_REGISTRY;

	/* read in object class */
	in >> surfclass;
	std::auto_ptr<Implicit> temp = NEW_IMPLICIT(surfclass);

	Implicit *imp = temp.release();
	if (!imp) {
		std::cerr << "Cannot find Surface of class " << surfclass;
		return in;
	}

	/* read in surface name */
	in >> c;
	if (c != '"') {
		std::cerr << "Specified object with class " << surfclass << "does not have a quoted name.\n";
		return in;
	}

	in >> std::noskipws;
	while ((in >> c) && c != '"') surfname += c;
	in >> std::skipws;

	imp->setObjectName(surfname);

	if (!(in >> c) || c != '{') {
		std::cerr << "Missing opening brace for surface " << surfname << " of type " << surfclass << ".\n";
		return in;
	}

	/* read general surface params */
	in >> &(imp->params);
	/* read specific surface params, like q's for implicit */
	imp->input(in);

	if (!(in >> c) || c != '}') {
		std::cerr << "Missing closing brace for surface " << surfname << " of type " << surfclass << ".\n";
		return in;
	}

	/* add into surfaces */
	surfaces->push_back(imp);

	/* eat up remaining whitespace */
	if (in >> c)
		in.putback(c);

	return in;
}


/// factory variables
/** Grab the Algebraics. */
extern int grabCone, grabCylinder, grabEllipsoid, grabTorus;

/** Grab the Geometrics. */
extern int grabPlane, grabPoint, grabSegment, grabSphere;

/** Grab the Operators. */
extern int grabOffset, grabImpList, grabSum, grabDifference, grabProduct;
extern int grabUnion, grabIntersection, grabMover, grabComplement;
extern int grabAUnion, grabAIntersection, grabASubtraction, grabCSGUnion, grabCSGIntersection, grabGradFdotV, grabSuggestiveSurface, grabZeroGaussianCurvature, grabRotatedImplicit;
extern int grabSpecularSurface; 
extern int grabDistanceField,grabKDOp;
extern int grabDirectionalDerivative,grabParabolicPoints,grabRadialCurvature;

#ifdef WB_USE_ITK
extern int grabITKImplicit;
#endif

#ifdef WB_USE_VTK
extern int grabVTKImplicit;
#endif

extern int grabConvolution, grabBlobby;

/** Grab the Blobs. */
extern int grabBlinn, grabWyvill;

/** Grab ADF and other wierd stuff. */
//extern int grabADF;
extern int grabRBF, grabInterfaceRBF, grabCompactRBF;
extern int grabSurfletCompactRBF, grabMultiscaleCSRBF, grabAdaptiveCSRBF, grabAdaptiveCSRBF_Precomp;

#ifdef WB_USE_SFL   
extern int grabAlexaImplicit, grabKolluriImplicit, grabProjNormImplicit;
#endif

extern int grabDummy;

extern int grabIcosahedral, grabThinPlateSpline;
extern int grabThinPlateSpline;

/** The compilation of this function activates implicit subclasses
 * in the factory.
 */
void grabby()
{
  grabCone++;
  grabCylinder++;
  grabEllipsoid++;
  grabTorus++;
  
  grabPlane++;
  grabPoint++;
  grabSegment++;
  grabSphere++;
  
  grabOffset++;
  grabImpList++;
  grabSum++;
  grabProduct++;
  grabDifference++;
  grabUnion++;
  grabIntersection++;
  grabMover++;
  grabComplement++;
  grabAUnion++;
  grabAIntersection++;
  grabASubtraction++;
  grabCSGUnion++;
  grabCSGIntersection++;
  
  grabGradFdotV++;
  grabSpecularSurface++;
  grabDirectionalDerivative++;
  grabZeroGaussianCurvature++;
  grabParabolicPoints++;
  grabRadialCurvature++;
  grabRotatedImplicit++;
#ifdef WB_USE_ITK  
  grabITKImplicit++;
#endif 
#ifdef WB_USE_VTK  
  grabVTKImplicit++;
#endif 
  grabSuggestiveSurface++;
  grabDistanceField++;
  grabKDOp++;
  
  grabConvolution++;
  grabBlobby++;

  grabBlinn++;
  grabWyvill++;

//  grabADF++;
  grabRBF++;
  grabInterfaceRBF++;
  grabCompactRBF++;
  grabSurfletCompactRBF++;
  grabMultiscaleCSRBF++;
  grabAdaptiveCSRBF++;
  grabAdaptiveCSRBF_Precomp++;
#ifdef WB_USE_SFL   
  grabAlexaImplicit++;
  grabKolluriImplicit++;
  grabProjNormImplicit++;
#endif

  grabIcosahedral++;
  grabThinPlateSpline++;
}
