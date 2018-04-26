/**
 * Declaration of a specular contour shader
 * @file ShaderSpecularContour.h
 * @date 04/07/2006
 * @author Matei N. Stroila
 * @remarks
 */

#ifndef SHADERSPECULARCONTOUR_H
#define SHADERSPECULARCONTOUR_H

#include "ShaderContour.h"

class LightPosition;

class ShaderSpecularContour : public ShaderContour
{
protected:	
	virtual void findTangent(const unsigned int i);	
public:
	
	MAKE_PARTICLESTUFF_NAME();
	
	// default constructor
	ShaderSpecularContour(Particles *ps=NULL);
	
	~ShaderSpecularContour();
	
private:

	LightPosition *light;
	double shine;
};

#endif