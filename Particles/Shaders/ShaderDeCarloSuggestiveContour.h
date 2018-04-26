/**
 * Declaration of ShaderDeCarloSuggestiveContour
 * @file ShaderDeCarloSuggestiveContour.h
 * @author Matei N. Stroila
 * @date  30/12/06
 */

#ifndef ShaderDeCarloSuggestiveContour_h
#define ShaderDeCarloSuggestiveContour_h

#include "ShaderContour.h"

class ParticleVelocity;
class LightPosition;

class ShaderDeCarloSuggestiveContour : public ShaderContour
{
	
protected:
	virtual void findTangent(const unsigned int i);
	virtual void checkCriticalPoint(){	};
	virtual void drawPost(void);
	
public:
		MAKE_PARTICLESTUFF_NAME();
	
	// default constructor
	ShaderDeCarloSuggestiveContour(Particles *ps=NULL);
	
	~ShaderDeCarloSuggestiveContour(){	};
	
private:
	
	LightPosition *light;
	bool useRadCurvNumerator;
	double dirDerivThreshold;
	
	double computeDerivative(const int i); //compute the directional derivative of the radial curvature in the view direction, at the position of the i-th particle

	ParticleVelocity *velocity;
	double velocityThreshold;
};

#endif
