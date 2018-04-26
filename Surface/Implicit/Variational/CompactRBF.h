/** @file CompactRBF.h
 * Declaration of an implicit surface class based on radial basis function interpolation.
 * Based off of the RBF.h file by John C. Hard and William Nagel.
 * @author: Mike Flavin
 * @date: April 13, 2004
 */
//Modified by Scott Kircher to make it work with surflet-style CSRBFs
//In particular, the KD-Tree has been completely replaced by a more
//generic (and hopefully more efficient) implementation.

#ifndef COMPACTRBF_H
#define COMPACTRBF_H

#include "Surface/Implicit/Variational/RBF.h"
#include "Surface/KDTree.h"


/** An implicit interpolation surface class based on radial basis functions.
 * Interpolating implicit surfaces can be generated from a collection
 * of radial basis functions. The resulting function is the sum
 * of these radial basis functions and will yield the assigned
 * value at the assigned constraint point. Hence the surface
 * interpolates the constraint values.
 *
 * Constraining points to zero values causes the surface to
 * interpolate the point. Constraining points to positive or
 * negative values defines interior/exterior and influences surface normals.
 *
 * This class implements the Wendland compactly-supported RBF for several different
 * desired levels of continuity for the implicit surface to be generated.
 *
 * Compactly supported RBFs have two main advantages - they make modeling with 
 * interpolating surfaces easier, since moving a control point has only a local
 * effect, and they also make the RBF matrix a sparse matrix, allowing the use
 * of fast sparse-matrix solution algorithms.
 */

class CompactRBF : public RBF
{
  public:

	/** Desired continuity for surface:
	 * Possible values:
	 * 0 = C^0
	 * 2 = C^2
	 * 4 = C^4
	 * 6 = C^6
	 * All other values are invalid, and will default to the nearest
	 * (lower) continuity (for example, 3 = C^2, 25 = C^6)
	 */
	int continuity;

	/** Support radius for RBFs.  Defaults to 1.
	 */
	double support_radius;

	/** Have the parameters of the interpolating surface changed since
	 *  we last solved the matrix?
	 */
	bool param_changed;

    /// Constructor.
    CompactRBF(void);
    
	/// Solve the constraint system to get the new interpolation surface.
    virtual void updateRBF(void);

	/// Functions for dealing with the parameter vector.
	unsigned int qlen();
	void getq(double *q);
	void _setq(double *q);
	void getqname(char **qn);

#ifndef INTERVAL_EVAL_ONLY
    virtual double proc(const gmVector3 & x);
    virtual gmVector3 grad(const gmVector3 & x);
    virtual gmMatrix3 hess(const gmVector3 & x);
#endif

	/// Functions for evaluating the surface and its derivatives.
    inline double phi(const gmVector3 & x);
    inline gmVector3 gradphi(const gmVector3 & x);
    inline gmMatrix3 hessphi(const gmVector3 & x);

/*Removed because I'm not sure they're correct, and we don't really need them --SIK
	inline Intervald phi(Box<double>);
    inline Box3d gradphi(Box<double>);
    inline IMatrix3d hessphi(Box<double>);
*/

	/// Function to generate a new interpolating surface from a set of constraints.
	virtual void interpolate(std::vector<gmVector3> positions, 
							 std::vector<double> interp_vals,
							 std::vector<gmVector3> normals,
							 std::valarray<bool>& flexible, 
							 bool pos_changed);

    MAKE_NAME();

protected:
	//KDtree for finding what points are close enough to have an influence
	KDTree<3,gmVector3> support_tree;

	/** Function to convert the modeler constraints to true RBF constraints
	 *   (in which each normal constraint is split into two separate RBF
	 *	  constraints, which contain simply a position and a desired function
	 *    value.)
	 */
	void convert_constraints();
	
	/// Regular versions of phi.
	inline double phi_c0(const gmVector3 & x);
	inline double phi_c2(const gmVector3 & x);
	inline double phi_c4(const gmVector3 & x);
	inline double phi_c6(const gmVector3 & x);

/*Removed because I'm not sure they're correct, and we don't really need them --SIK
	/// Interval versions of phi.
	inline Intervald phi_c0(Box<double> b);
	inline Intervald phi_c2(Box<double> b);
	inline Intervald phi_c4(Box<double> b);
	inline Intervald phi_c6(Box<double> b);
*/

	/// Regular versions of gradphi.
	inline gmVector3 gradphi_c0(const gmVector3 & x);
	inline gmVector3 gradphi_c2(const gmVector3 & x);
	inline gmVector3 gradphi_c4(const gmVector3 & x);
	inline gmVector3 gradphi_c6(const gmVector3 & x);

/*Removed because I'm not sure they're correct, and we don't really need them --SIK
	/// Interval versions of gradphi.
	inline Box3d gradphi_c0(Box<double>b);
	inline Box3d gradphi_c2(Box<double>b);
	inline Box3d gradphi_c4(Box<double>b);
	inline Box3d gradphi_c6(Box<double>b);
*/

	/// Regular versions of hessphi.
	inline gmMatrix3 hessphi_c0(const gmVector3 & x);
	inline gmMatrix3 hessphi_c2(const gmVector3 & x);
	inline gmMatrix3 hessphi_c4(const gmVector3 & x);
	inline gmMatrix3 hessphi_c6(const gmVector3 & x);

/*Removed because I'm not sure they're correct, and we don't really need them --SIK
	/// Interval versions of hessphi.
	inline IMatrix3d hessphi_c0(Box<double> b);
	inline IMatrix3d hessphi_c2(Box<double> b);
	inline IMatrix3d hessphi_c4(Box<double> b);
	inline IMatrix3d hessphi_c6(Box<double> b);
*/

};

#endif

