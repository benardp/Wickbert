/**
 * Declaration of a shader for silhouette contours
 * @file ShaderSilhouetteContour.h
 * @date 6/25/2005
 * @author Matei N. Stroila
 * @remarks
 */
#ifndef SHADERSILHOUETTESIMPLE_H
#define SHADERSILHOUETTESIMPLE_H


#include "ShaderContourSimple.h"

class ShaderSilhouetteSimple : public ShaderContourSimple
{	
protected:	

	/** Compute the tangent to the contour at the location of the i-th particle. 
	 * @param i the particle index
	 */
	virtual void findTangent(const unsigned int i);
	virtual void drawPost();
	
public:
		
	MAKE_PARTICLESTUFF_NAME();
	
	// default constructor
	ShaderSilhouetteSimple(Particles *ps=NULL);
	
	~ShaderSilhouetteSimple(){};
	
	void attachAttributes();
};

#endif

