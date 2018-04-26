/**
 * Declaration of a shadow contour shader 
 * @file ShaderShadowContour.h
 * @date 06/25/2005
 * @author Matei N. Stroila
 * @remarks
 */

#ifndef SHADERSHADOWCONTOUR_H
#define SHADERSHADOWCONTOUR_H

#include "ShaderContour.h"

class LightPosition;

class ShaderShadowContour : public ShaderContour
{
protected:	
	virtual void findTangent(const unsigned int i);
	virtual void checkCriticalPoint4D(double *pstatus);
	virtual void checkCriticalPoint5D(double *pstatus);
	virtual void checkCriticalPointPC(); //check for critical points using the parabolic contours 
	
public:
	
	MAKE_PARTICLESTUFF_NAME();
	
	// default constructor
	ShaderShadowContour(Particles *ps=NULL);
	~ShaderShadowContour();

	void getShadeColor(std::vector<short>& shadeColor) const{
		shadeColor.push_back((short) (shadeColorf[0] * 255.));
		shadeColor.push_back((short) (shadeColorf[1] * 255.));
		shadeColor.push_back((short) (shadeColorf[2] * 255.));
	}
	
private:

	LightPosition *light;
	gmVector3 shadeColorf;
	struct CPparams p;
	short shadeColor[3];

};

#endif