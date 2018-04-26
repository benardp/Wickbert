/** @file AdaptiveCSRBF.h
 * Extension of SurfletCompactRBF to adaptively choose CSRBF centers
 * rather than have them be at every input data point. Following Ohtake, et al...
 * "3D Scattered Data Approximation with Adaptive Compactly Supported Radial Basis Functions"
 * Only the Adaptive Partition of Unity is implemented, since the RBF fitting part
 * would have almost no impact on the particle visualization (and the adaptive PU is slow enough already).
 * @author: Scott Kircher
 * @date: November 28, 2004
 */

#ifndef ADAPTIVECSRBF_H
#define ADAPTIVECSRBF_H

#include "Surface/Implicit/Variational/SurfletCompactRBF.h"

/** AdaptiveCSRBF class
 */
class AdaptiveCSRBF : public SurfletCompactRBF
{
public:
	MAKE_NAME();

	//default constructor
    AdaptiveCSRBF();

	/// Solve the constraint system to get the new interpolation surface.
    virtual void updateRBF(void);

	unsigned int qlen();
	void getq(double *q);
	void _setq(double *q);
	void getqname(char **qn);

#ifndef INTERVAL_EVAL_ONLY
    virtual double proc(const gmVector3 & x);
    virtual gmVector3 grad(const gmVector3 & x);
    virtual gmMatrix3 hess(const gmVector3 & x);
#endif

	/** @name AdaptivePhi
	 * Functions for evaluating the surface and its derivatives.
	 * Unlike base classes SurfletCompactRBF and CompactRBF, we
	 * allow a different support radius at every CSRBF center.
	 * @{ */
    inline double phi(const gmVector3 & x,double support);
    inline gmVector3 gradphi(const gmVector3 & x,double support);
    inline gmMatrix3 hessphi(const gmVector3 & x,double support);
	/// @}

protected:
	//sparsity parameter (the higher the value, the larger the support sizes will be, and the 
	//fewer centers will be used)
	double sparsity;

	//in addition to the centers and surflets we already have (from CompactRBF and SurfletCompactRBF respectively)
	//we also need a support size and a density weight at every center location
	std::vector<double> local_support_size;
	std::vector<double> local_density_weight;

	//diagonal length (not necessarily main diagonal... we assume data's BB is roughly cube-shaped
	double diagonal;

	//overlap threshold (how much each point should be "covered" by the approximations
	double overlap_threshold;

	//the minimum allowed support radius
	double min_support_size;
	
	//maximum support radius used (size of radius query for KD-tree)
	double max_support_size; //not adjustable, this is measured

	//once the approximation centers have been added, we don't want to have them all in our KDtree
	//(we want only the active ones). So, we need a map from support-tree indices to constraint indices
	std::vector<int> approximation_indices;

	/** Function to convert the modeler constraints to true Surflet RBF constraints
	 * unlike RBF and CompactRBF, SurfletCompactRBFs do not add extra centers for
	 * normal constraints, instead they fit a surflet to the surface at that point.
	 */
	void convert_constraints(); //this includes building the KDtree

	// For finding optimal local support sizes
	double computeLocalWeights(const gmVector3 &point,double sigma,std::vector<int> &nz_centers,std::vector<double> &nz_weights);
	double sparseApproxError(const gmVector3 &point,const gmVector3 &normal,double sigma);
	double findOptimalSupportSize(const gmVector3 &point,const gmVector3 &normal);
};

#endif