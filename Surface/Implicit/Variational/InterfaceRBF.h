/** @file RBF.h
 * Declaration of an implicit surface class based on radial basis function interpolation.
 * It offers the possibility to exchange the function phi (important for smoothing Carr et al. 2003). 
 * Also calculation of the RBF can done based on particles
 * Therefore it is possible to create a particles object whose position and values can be connected.
 * E.g.: Create a particle system that spreads out over the surface based on curvature and then use these 
 * particles to calculate an appropriate RBF. Or place particles on a mesh etc.
 * @author: Elmar Eisemann and I some parts from RBF to not duplicate code. ;)
 * @date: Sommer, 2006
**/

#ifndef INTERFACERBF_H
#define INTERFACERBF_H

#include "Surface/SurfParam.h"
#include "Surface/Implicit/Variational/RBF.h"


//class Particles;
class ParticlePosition;
class ParticleScalar;
class RBFBasicFunction;



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
class InterfaceRBF : public RBF
{
	//callback to interpolate the RBF from the given particle values
	class InterpolateParticlesWithRBF : public SurfParamButton::Callback
	{
	public:
		InterfaceRBF *pbRBF;
		InterpolateParticlesWithRBF(InterfaceRBF * pp) {pbRBF = pp;}
		virtual void onbuttonpress() {pbRBF->interpolate();}
	};

	//callback to load the RBF 
	class LoadRBFCallback : public SurfParamButton::Callback
	{
	public:
		InterfaceRBF *pbRBF;
		LoadRBFCallback(InterfaceRBF * pp) {pbRBF = pp;}
		virtual void onbuttonpress() {pbRBF->load();}
	};

	//callback to save the RBF
	class SaveRBFCallback : public SurfParamButton::Callback
	{
	public:
		InterfaceRBF *pbRBF;
		SaveRBFCallback(InterfaceRBF * pp) {pbRBF = pp;}
		virtual void onbuttonpress() {pbRBF->save();}
	};

	ParticlePosition * zeroConstraints;

	//this function is not initialized, nor allocated here
	//it is given as a paramter to the RBFBasicFunction selector
	//that class keeps track of this allocation
	//It would have been more beautiful to have only a pointer to
	//the selector and query the phiFunction but this would 
	//have implied one supplementary redirection.
	RBFBasicFunction * _phiFunction;
	std::string _filename;
public:
    /// Constructor.
    InterfaceRBF(void);

	void load();
	void save() const;

	/// Overriden to read in constraints
	/// we should possibly one day also save the type of function that was used
	virtual bool readImplicit(std::ifstream &file, bool verbose);

	virtual void interpolate();
	virtual void updateRBF();

    virtual double proc(const gmVector3 & x);
    virtual gmVector3 grad(const gmVector3 & x);
    virtual gmMatrix3 hess(const gmVector3 & x);
	
	/// Functions for dealing with the parameter vector, here it is actually empty
	virtual unsigned int qlen();
	virtual void getq(double *q);
	virtual void getqname(char **qn);

	virtual void _setq(double *q);
	//there is no parameter derivative
	void procq(const gmVector3 & x,double *dq);


	virtual Intervald proc(const Box<double>& );
    virtual Box3d grad(const Box<double>& );
    virtual IMatrix3d hess(const Box<double>& );

    MAKE_NAME();
};

#endif 

