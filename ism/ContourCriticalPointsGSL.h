/**
 * Declaration of a class used to find the critical points of a feature contour
 * using GSL solvers.
 * @file ContourCriticalPointsGSL.h
 * @date 09/05/2005
 * @author Matei N. Stroila
 * @remarks
 */
#ifndef CONTOURCRITICALPOINTSGSL_H
#define CONTOURCRITICALPOINTSGSL_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include "Surface/Box.h"
#include "Newton.h"

class Implicit;
class RotatedImplicit;
class SingularityClassification;


struct CPparams
{
	//pointer to implicit
    Implicit* F;
	//camera zoom
    double radius;
	//light position (if needed)
	double Lx;
	double Ly;
	double Lz;
	double phi0, phi1, theta0, theta1; //Euler angles at t0 and t1
	RotatedImplicit* pRF;
	gmMatrix3* pRinit; //pointer to initial rotation
	gmMatrix3* pRfinal; //pointer to final roatation
};


static int epsilon[3][3][3] = {{{0,0,0},{0,0,1},{0,-1,0}},
	{{0,0,-1},{0,0,0},{1,0,0}},
	{{0,1,0},{-1,0,0},{0,0,0}}};


typedef enum{ SILHOUETTES=1, SHADOWS}
ContourType;

class ContourCriticalPointsGSL{
	
public: 
	
	ContourCriticalPointsGSL(CPparams *params, double _NewtonTolerance);
	
	~ContourCriticalPointsGSL();

	/**
	 * Function that returns true if no zeros are in the Box X and false if cannot decide it.
	 * @param X box of intervals where the function checks for possible solutions.
	 */
	bool NewtonSubdivisionBreak4D(const Box<double> & X, void *params, ContourType type);
	bool NewtonSubdivisionBreak5D(const Box<double> & X, void *params, ContourType type);
	bool NewtonSubdivisionBreak_kGkR_4D(const Box<double> & X, void *params, ContourType type);
	
	/**
	 * critical point finding using interval subdivision to eliminate large areas
	 & and GSL solvers to find the zeros in smaller boxes.
	 */
	int FindZeros4D(Box<double> &X, void *params, ContourType type, ClassifiedPointList& cpList, double* status);
	int FindZeros5D(Box<double> &X, void *params, ContourType type, ClassifiedPointList& cpList, double* pstatus);
	int FindZeros_kGkR_4D(Box<double> &X, void *params, ContourType type, ClassifiedPointList& cpList, double* pstatus);
	

	/**
	 * Interval arithmetic based critical point finding.
	 */
	int FindZerosIA(Box<double> &X, void *params, ContourType type, ClassifiedPointList& cpList, double* pstatus);	
	int FindZeros5IA(Box<double> &X, void *params, ContourType type, ClassifiedPointList& cpList, double* pstatus);
	
	/**
	 * Solvers using GSL.
	 */
	int solveContourCriticalPointGSL5D( Box<double>& theBoxBounds, void *params, ContourType type, ClassifiedPointList& cpList);
	int solveContourCriticalPointGSL4D(Box<double>& theBoxBounds, void *params, ContourType type, ClassifiedPointList& cpList);
	int solveContourCriticalPointGSL_kGkR_4D(Box<double>& theBoxBounds, void *params, ContourType type, ClassifiedPointList& cpList);
	
	/**
	 * Function that returns true if the root newroot is already in the root list cpList.
	 * @param newroot New critical point.
	 * @param cpList List of critical points already found.
	 */
	bool SameRoot(TNTPoint& newroot, ClassifiedPointList& cpList);
	
	/**
	 * Helper function that outputs the state of the GSL solver
	 */
	void print_state (size_t iter, gsl_multiroot_fdfsolver * s);
	
	/**
	 * The maximum width of a box that is no longer subdivided. A solver is then used to find the roots.
	 */
	double NewtonTolerance;
	/**
	 * The maximum width of a box that is considered to be a root.
	 */
	double WidthTolerance;
	/**
	 * Threshold value for the distance between two critical points to be considered the same critical point.
	 */
	double SameTolerance;
	/**
	 * The maximum number of iterations in a GSL solver.
	 */
	unsigned int maxIter;
	
	/**
	 * Function that computes the 4D system of 
	 * Simon Plantinga, Gert Vegter "Contour generators of evolving implicit surfaces" 
	 * http://doi.acm.org/10.1145/781606.781614 for shadows.
	 * @param x the 4D point.
	 * @param f the 4D value of the function at x.
	 */
	static int shadowContour4DCP_f (const gsl_vector * x, void *params, gsl_vector * f);
	/**
	 * Function that computes the Jacobian of the 4D system of 
	 * Simon Plantinga, Gert Vegter for shadows.
	 * @param x the 4D point.
	 * @param J the Jacobian at x.
	 */
	static int shadowContour4DCP_df (const gsl_vector * x, void *params,  gsl_matrix * J);
	/**
	 * Function that puts together (for GSL) f and the Jacobian of the 4D system of 
	 * Simon Plantinga, Gert Vegter for shadows.
	 * @param x the 4D point.
	 * @param f the 4D value of the function at x.
	 * @param J the Jacobian at x.
	 */
	static int shadowContour4DCP_fdf (const gsl_vector * x, void *params, gsl_vector * f, gsl_matrix * J);
	/**
	 * Interval version of @see shadowContour4DCP_f().
	 * @param X the 4D box.
	 * @param F the 4D box of the function at X.
	 */
	static int shadowContour4DCP_F(const Box<double> & X, void *params, Box<double>& F);
	/**
	 * Function that computes the 4D system of 
	 * Simon Plantinga, Gert Vegter for silhouettes.
	 * @param x the 4D point.
	 * @param f the 4D value of the function at x.
	 */
	static int
		silhouetteContour4DCP_f (const gsl_vector * x, void *params, 
								 gsl_vector * f);
	/**
	 * Interval version of @see silhouetteContour4DCP_f().
	 * @param X the 4D box.
	 * @param Fi the i-th interval of the function at X.
	 */
	void silhouetteContour4DCP_F(const Box<double> & X, void *params, Interval<double>& Fi, const int i);
	

	/**
	 * Function that computes the 5D system for silhouettes.
	 * @param x the 5D point.
	 * @param f the 5D value of the function at x.
	 */
	static int contourCP_f (const gsl_vector * x, void *params, 
					 gsl_vector * f);
	/**
	 * Interval version of @see contourCP_f()
	 * @param X the 5D box.
	 * @param F the 5D box of the function at X.
	 */
	static int contourCP_F(const Box<double> & X, void *params, Box<double>& F);
	/**
	 * Function that computes the Jacobian of the 5D system for silhouettes.
	 * @param x the 5D point.
	 * @param J the Jacobian at x.
	 */
	static int contourCP_df (const gsl_vector * x, void *params, 
					  gsl_matrix * J);
	/**
	 * Interval version of @see contourCP_df()
	 * @param X the 5D box.
	 * @param Jacobian the 5D interval matrix at X.
	 */
	static int contourCP_DF(const Box<double> & X, void *params, IMatrix& J);
	
	/**
	 * Function that puts together (for GSL) f and the Jacobian of the 5D system 
	 * for silhouettes.
	 * @param x the 5D point.
	 * @param f the 5D value of the function at x.
	 * @param J the Jacobian at x.
	 */
	static int contourCP_fdf (const gsl_vector * x, void *params,
					   gsl_vector * f, gsl_matrix * J);
	/**
	 * Function that computes the 5D system for shadows.
	 * @param x the 5D point.
	 * @param f the 5D value of the function at x.
	 */
	static int shadowContourCP_f (const gsl_vector * x, void *params, gsl_vector * f);
	/**
	 * Interval version of @see shadowContourCP_f()
	 * @param X the 5D box.
	 * @param F the 5D box of the function at X.
	 */
	static int shadowContourCP_F(const Box<double> & X, void *params, Box<double>& F);
	/**
	 * Interval version of @see shadowContourCP_df()
	 * @param X the 5D box.
	 * @param Jacobian the 5D interval matrix at X.
	 */
	static int shadowContourCP_DF(const Box<double> & X, void *params, IMatrix& Jacobian);
	/**
	 * Function that computes the Jacobian of the 5D system for shadows.
	 * @param x the 5D point.
	 * @param J the Jacobian at x.
	 */
	static int shadowContourCP_df (const gsl_vector * x, void *params, gsl_matrix * J);
	/**
	 * Function that puts together (for GSL) f and the Jacobian of the 5D system 
	 * for shadows.
	 * @param x the 5D point.
	 * @param f the 5D value of the function at x.
	 * @param J the Jacobian at x.
	 */
	static int shadowContourCP_fdf (const gsl_vector * x, void *params, gsl_vector * f, gsl_matrix * J);
	
	/**
	 * Function that computes the 4D system (F = Fv = kG = kR = 0) for silhouettes.
	 * @param x the 4D point.
	 * @param f the 4D value of the function at x.
	 */
	static int silh_kGkR_4DCP_f (const gsl_vector *x, void *params, gsl_vector *f);
	/**
	 * Function that computes the Jacobian of the 4D system (F = Fv = kG = kR = 0) 
	 * for silhouettes.
	 * @param x the 4D point.
	 * @param J the Jacobian at x.
	 */
	static int silh_kGkR_4DCP_df (const gsl_vector * x, void *params, gsl_matrix * J);
	/**
	 * Function that puts together (for GSL) f and the Jacobian of the 4D system 
	 * F = Fv = kG = kR = 0) for silhouettes.
	 * @param x the 4D point.
	 * @param f the 4D value of the function at x.
	 * @param J the Jacobian at x.
	 */
	static int silh_kGkR_4DCP_fdf (const gsl_vector * x, void *params, gsl_vector * f, gsl_matrix * J);
	/**
	 * Interval version of @see silh_kGkR_4DCP_f(). Compute only the ith component of the system.
	 * @param X the 4D box.
	 * @param Fi the interval of the ith function in the system at X.
	 */
	void silh_kGkR_4DCP_F(const Box<double>& X, void *params, Interval<double>& Fi, const int i);	
	/**
	 * Function that computes the 4D system (F = Fl = kG = kR = 0) for shadows.
	 * @param x the 4D point.
	 * @param f the 4D value of the function at x.
	 */
	static int shadow_kGkR_4DCP_f (const gsl_vector *x, void *params, gsl_vector *f);
	/**
	 * Interval version of @see shadow_kGkR_4DCP_f().
	 * @param X the 4D box.
	 * @param F the 4D box of the function at X.
	 */
	static int shadow_kGkR_4DCP_F(const Box<double>& X, void *params, Box<double>& F);
	
	
private:

	Implicit *myF;
	gmMatrix3 Rinit, Rfinal; //initial (t0) and final (t1) rotation matrices
	void  computeInitialFinalRotations(void *params);
	RotatedImplicit* pRF;
	SingularityClassification* singClassif;
	double customWidth(const Box<double> & X);

};	
	

#endif
