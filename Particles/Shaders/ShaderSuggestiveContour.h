/**
 * Declaration of ShaderSuggestiveContour
 * @file ShaderSuggestiveContour.h
 * @author Matei N. Stroila
 * @date  7/16/05
 */

#ifndef SHADERSUGGESTIVECONTOUR_H
#define SHADERSUGGESTIVECONTOUR_H

#include "ShaderContour.h"

class ParticleVelocity;
class LightPosition;

class ShaderSuggestiveContour : public ShaderContour
{
	
protected:
	virtual void findTangent(const unsigned int i);
	virtual void checkCriticalPoint(){	};
	virtual void drawPost(void);
	
public:
		MAKE_PARTICLESTUFF_NAME();
	
	// default constructor
	ShaderSuggestiveContour(Particles *ps=NULL);
	
	~ShaderSuggestiveContour(){	};
	
private:
	
	LightPosition *light;
	
	double dirDerivThreshold;
	
	double computeDerivative(const int i); //compute the directional derivative of the radial curvature in the view direction, at the position of the i-th particle

	ParticleVelocity *velocity;
	double velocityThreshold;
};

#endif
