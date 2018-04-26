/** @file RBF.h
 * Declaration of an implicit surface class based on radial basis function interpolation.
 * @author: John C. Hart, William Nagel
 * @date: Fall Semester, 2000, updated 23 Sep 2003
 */

#ifndef RBF_H
#define RBF_H

#include "Surface/Implicit/Implicit.h"

//class Particles;
class ParticlePosition;
class ParticleScalar;

class RBFPoint 
{
public:
    /// Location of constraint point
    gmVector3 c;

    /// Desired function value at point c.
    double h;

    /** Weight of constraint point.
     * Initially unknown, and solved for later.
     */
    double w;
	
	RBFPoint(gmVector3 c, double h) { this->c = c; this->h = h; }

	RBFPoint() {c = gmVector3();}
};

class RBFModelerConstraint
{
public:
	/// Location of constraint
	gmVector3 c;

	/// Desired function value at point c.
	double h;

	/// Desired gradient (normal) at point c.
	gmVector3 n;

	/// Is the current constraint a normal constraint or not?
	bool nc;

	RBFModelerConstraint(gmVector3 c,  gmVector3 n)
	{
		n.normalize();

		this->c = c;
		this->n = n;
		this->h = 0.0;
		this->nc = true;
	}

	RBFModelerConstraint(gmVector3 c, double h)
	{
		this->c = c;
		this->n = gmVector3(0.0, 0.0, 0.0);
		this->h = h;
		this->nc = false;
	}

	RBFModelerConstraint(gmVector3 c, const gmVector3& n, double h, bool nc)
	{
		this->c = c;
		this->n = n;
		this->h = h;
		this->nc = nc;
	}

	RBFModelerConstraint() {c = gmVector3();}
};

typedef std::vector<RBFPoint> RBFPoints;
typedef std::vector<RBFModelerConstraint> RBFModelerConstraints;

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
 */
class RBF : public Implicit
{
  public:
#if 0
	/** Pointer to a Particles object in the wickbert modeler.  Intended so that
	* an RBF object can take in a Particles object as its list of centers 
	*/
	Particles *p;
	std::string filename;
#else
	ParticlePosition *positions;
	ParticleScalar *values;
	ParticleScalar *weights;
#endif

	/** Pointer to an RBFPoints vector that constrain the implicit.
	 * Implemented as a pointer so different RBF classes and subclasses
	 * can share the same data via the RBFPoints class interface.
	 */
	RBFPoints centers;
	
	/** An RBFModelerConstraints vector that contains the "conceptual" view
	 * of the constraints (in other words, a normal constraint is a constraint
	 * with a position and desired normal, not two nearby value constraints).
	 */
	RBFModelerConstraints constraints;
	  
	/** Affine part of radial basis function surface.
     */
    double m_p[4];

    /// Constructor.
    RBF(void);

	/// Overriden to read in constraints
	virtual bool readImplicit(std::ifstream &file, bool verbose);

    /// Solve the constraint system to get the new interpolation surface.
    virtual void updateRBF(void);

	//we don't want setq being called all the time
	//to do finite differences on q... and furthermore
	//gradient of q is a little meaningless for these kinds fo surfaces
	//--SIK
	void procq(const gmVector3 & x,double *dq);
//	{
//		for(int c=0;c<qlen();c++) dq[c]=0.0;
//	}

    
#ifndef INTERVAL_EVAL_ONLY
    virtual double proc(const gmVector3 & x);
    virtual gmVector3 grad(const gmVector3 & x);
    virtual gmMatrix3 hess(const gmVector3 & x);
#endif
	
	/// Functions for dealing with the parameter vector.
	virtual unsigned int qlen();
	virtual void getq(double *q);
	virtual void getqname(char **qn);

	virtual void _setq(double *q);
    
	virtual Intervald proc(const Box<double>& );
    virtual Box3d grad(const Box<double>& );
    virtual IMatrix3d hess(const Box<double>& );

    virtual double phi(const gmVector3 & x);
    virtual gmVector3 gradphi(const gmVector3 & x);
    virtual gmMatrix3 hessphi(const gmVector3 & x);

    virtual Intervald phi(const Box<double>& );
    virtual Box3d gradphi(const Box<double>& );
    virtual IMatrix3d hessphi(const Box<double>& );

//	virtual void interpolate(std::vector<gmVector3> positions, std::vector<double> interp_vals, std::vector<gmVector3> normals, std::valarray<bool>& flexible, bool pos_changed);
	//virtual void interpolate(Particles *ps, float phi);
	virtual void interpolate();

    MAKE_NAME();

private:
	/** Have the parameters of the interpolating surface changed since
	 *  we last solved the matrix?
	 */
	bool param_changed;
};

#endif 

