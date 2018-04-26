/**
 * Declaration of a shader for silhouette contours
 * @file ShaderIntersectionContour.h
 * @date 6/25/2005
 * @author Matei N. Stroila
 * @remarks
 */
#ifndef SHADERINTERSECTIONCONTOUR_H
#define SHADERINTERSECTIONCONTOUR_H


#include "ShaderContourSimple.h"

class ShaderIntersectionContour : public ShaderContourSimple
{	
protected:	

	/** Compute the tangent to the contour at the location of the i-th particle. 
	 * @param i the particle index
	 */
	virtual void findTangent(const unsigned int i);
	virtual void drawPost();
	
public:
		
	MAKE_PARTICLESTUFF_NAME();

	ImplicitInterrogator *secondsurf;

	// default constructor
	ShaderIntersectionContour(Particles *ps=NULL);
	
	~ShaderIntersectionContour(){};
	
	void attachAttributes();
};

#endif

