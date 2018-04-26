/**
 * Declaration of the particle visibility attribute.
 * @file ParticleVisibility.h
 * @author Matei N. Stroila
 * @date 6/10/2005
 * @remarks 
 */

#ifndef PARTICLEVISIBILITY_H
#define PARTICLEVISIBILITY_H

#include "ParticleAttribute.h"

class ParticlePosition;
class ParticleNormal;
class ViewDependence;
class ParticleBoundingBox;
class ImplicitInterrogator;
class Implicit;

/*! \class ParticleVisibility
*  \brief ParticleVisibility is a ParticleAttribute which represents the
*         visibility of particles in a system with a vector of doubles.
*/

class ParticleVisibility : public ParticleAttribute
{
	friend class VisibilityMethodCallback;

public:
    MAKE_PARTICLESTUFF_NAME();
	
	/// Constructor adds particle visibility to a system of particles.
	ParticleVisibility(Particles *ps=NULL, const std::string& name=std::string("ParticleVisibility"));
	
	~ParticleVisibility(){};
	virtual void prepare();

	void clear();
	void reset();
	
    /*! This method returns the visibility of the ith particle.
		*  \param i The particle of interest.
		*  \return mVisibility[i]
		*/
	
    inline const int getVisibility(const unsigned int i) const
    {
		return mVisibility[i];
    } // end ParticleVisibility::getVisibility()
	
    /*! OVERRIDDEN!!! THIS FUNCTION IS NO LONGER INFLUENCED BY THE CAMERA!
		*  \param i The particle of interest.
		*  \return mVisibility[i]
		*/
	int getVisibility( gmVector3* cameraPosition,  unsigned int i);
	
	/// Callback for particle addition.
	virtual void particleAdded();
	
	/// Callback for particle removal.
	virtual void particleRemoved(unsigned int i);
	   
	inline void setVisibility(const unsigned int i, const int j)
    {
		mVisibility[i] = j;
    } 
	
	inline void setOffset(double o)
    {
		offset = o;
    }
	
	inline void setUseNormalsTest(bool b)
    {
		useNormalTest = b;
    }


	
private:
	enum ParticleVisibilityMethod {ParticleVisibility_GSL, ParticleVisibility_STEPPING} visibilityMethod;

	class VisibilityMethodCallback : public PSParamComboBox::Callback
	{
		ParticleVisibility * _pv;
	public:
		VisibilityMethodCallback(ParticleVisibility * pv);
		void itemselected(unsigned int selection);
	};
	int computeVisibilityGSL(const int i) const;
	int computeVisibilityStepping(const int i) const;

	/// ParticleVisibility encodes yes/no visibility for each particle in a particle system.
	std::vector<int> mVisibility;
	std::vector<Implicit*> _implicits;

	ParticlePosition *pos;
	ViewDependence *viewDep;
	ParticleNormal *orient;	
	ParticleBoundingBox* bbox;
	gmVector3 oldCameraPosition;


	double newtonEps; //epsilon for Newton method
	int maxIter; //max number of Newton iterations
	double offset; //particle offset position to avoid tangential intersections (for silhouettes)
	double initX; //initial x for iterations
	double minX; //minimum x
	double bboxOffset, nbStepsPerUnit;
	bool useNormalTest, useBbox;

	bool lazyEvaluation;
	//don't evaluate visibility
	bool deactivated;
	//temp variable indicating whether we need to reinit
	bool isDeactivated;
  /**
	* Compute the intersection @param intersectBbox of the view vector with the bounding box
	*/
	void computeIntersection(const int i, gmVector3& intersectBbox);
	
	
}; // end class ParticleVisibility

//Global functions for GSL
double myF (double x, void *params);
double myF_deriv (double x, void *params);
void myF_fdf (double x, void *params, 
				 double *y, double *dy);
struct myF_params
{
    Implicit *myImp;
	gmVector3 *myIntersectBbox;
	gmVector3 *myParticlePos;
	ParticleBoundingBox* bbox;
};

#endif // PARTICLEVISIBILITY_H
