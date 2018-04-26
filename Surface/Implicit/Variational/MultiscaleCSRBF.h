/** @file MultiscaleCSRBF.h
 * Uses a hierarchy of SurfletCompactRBFs to interpolate a set of points
 * in a multiscale manner.
 * @author: Scott Kircher
 * @date: November 7, 2004
 */

#ifndef MULTISCALECSRBF_H
#define MULTISCALECSRBF_H

#include "Surface/Implicit/Variational/SurfletCompactRBF.h"

/** MultiscaleCSRBF class
  */
class MultiscaleCSRBF : public RBF
{
public:
    MAKE_NAME();

/***Set of parameters***/
protected:
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
	int max_level;
	float support_scale; //multiplier on support size
	//surflet degree
	//valid values:
	//  1 - Local linear approximations
	//  2 - Local quadratic approximations
	int surflet_degree;

public:
	/// Functions for dealing with the parameter vector.
	unsigned int qlen();
	void getq(double *q);
	void _setq(double *q);
	void getqname(char **qn);
/***End parameters***/

	//default constructor
    MultiscaleCSRBF();
	//destructor
	virtual ~MultiscaleCSRBF() {clearLevels();}

	void clearLevels();

	/// Solve the constraint system to get the new interpolation surface.
    virtual void updateRBF(void);

#ifndef INTERVAL_EVAL_ONLY
    virtual double proc(const gmVector3 & x);
    virtual gmVector3 grad(const gmVector3 & x);
    virtual gmMatrix3 hess(const gmVector3 & x);
#endif

	/// Function to generate a new interpolating surface from a set of constraints.
	virtual void interpolate(std::vector<gmVector3> positions, 
							 std::vector<double> interp_vals,
							 std::vector<gmVector3> normals,
							 std::valarray<bool>& flexible, 
							 bool pos_changed);

protected:
	//vector of CSRBF levels
	std::vector<SurfletCompactRBF *> levels;

	bool param_changed;
};

#endif

