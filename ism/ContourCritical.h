/**
 * Declaration of a class used to find the critical points of a feature contour
 * using pure interval arithmetic methods.
 * @file ContourCritical.h
 * @date 01/13/2006
 * @author Matei N. Stroila
 * @remarks
 */
#ifndef CONTOURCRITICAL_H
#define CONTOURCRITICAL_H

#include "Newton.h"
#include "ContourCriticalPointsGSL.h"

class ContourCritical : public virtual Newton
{
	
public:
    ContourCritical();
	~ContourCritical();
    ContourCritical(CPparams*, double _NewtonTolerance, bool _useNewton);	
	int search(Box<double> &X, ContourType type, ClassifiedPointList& cpList);
	
protected:
	/// Overridden NewtonEquation to solve A(Y-Xc) = b for critical points.
	virtual bool NewtonEquation(Box<double> &X, IMatrix &A, Box<double> &b);
    /// Overridden Newton solver method to stop subdividing.
    virtual bool NewtonSubdivisionBreak(Box<double> & X);
	
private:
		
	CPparams *params;
	/**
	 * Pointer to an implicit object that is the rotated implicit in the camera coordinate system.
	 * prF is used at the classification of the critical points. 
	 */
	RotatedImplicit *pRF;
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
	 * Initial (at t0) and final (at t1) rotation matrices.
	 */
	gmMatrix3 Rinit, Rfinal; 
	/**
	 * Pointer to an object used to classify the critical points
	 */
	SingularityClassification* singClassif;
	
	/**
	 * Function that computes the 4D system (F = Fv = kG = kR = 0) for silhouettes.
	 * @param x the 4D point.
	 * @param f the 4D value of the function at x.
	 */
	void silh_kGkR_4DCP_f (const TNTPoint& x, TNTPoint &f);
	/**
	 * Interval version of @see silh_kGkR_4DCP_f(). Compute only the ith component of the system.
	 * @param X the 4D box.
	 * @param Fi the interval of the ith function in the system at X.
	 */
	void silh_kGkR_4DCP_F(const Box<double>& X, Interval<double>& Fi, const int i);
	/**
	 * Function that computes the interval matrix Jacobian of the 4D system (F = Fv = kG = kR = 0) 
	 * for silhouettes.
	 * @param X the 4D box.
	 * @param J the Jacobian at X.
	 */
	void silh_kGkR_4DCP_DF(const Box<double> & X, IMatrix& J);
	void silh_interp_f (const TNTPoint& x, TNTPoint &f);
	void silh_interp_F(const Box<double>& X, Interval<double>& Fi, const int i);
	void silh_interp_DF(const Box<double> & X, IMatrix& Jacobian);
	
};

#endif

