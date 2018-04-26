/**
 * Declaration of a shader for silhouette contours
 * @file ShaderSilhouetteContour.h
 * @date 6/25/2005
 * @author Matei N. Stroila
 * @remarks
 */
#ifndef SHADERSILHOUETTECONTOUR_H
#define SHADERSILHOUETTECONTOUR_H


#include "ShaderContour.h"

class ShaderSilhouetteContour : public ShaderContour
{	
protected:	

	/** Compute the tangent to the contour at the location of the i-th particle. 
	 * @param i the particle index
	 */
	virtual void findTangent(const unsigned int i);
	/** Check for critical points using the Plantinga Vegter system. 
	 */
	virtual void checkCriticalPoint4D(double *pstatus);
	/** Check for critical points using a 5D system. 
	 */
	virtual void checkCriticalPoint5D(double *pstatus);
	/** Check for critical points using the parabolic contours. 
	 */
	virtual void checkCriticalPointPC(); 
	/** Check for critical points using the Gaussian and the radial curvatures. 
	 */
	virtual void checkCriticalPoint_kGkR_4D(double *pstatus);
	/** Check for critical points using a pure interval method. 
	 */
	virtual void checkCriticalPointInterval();

	virtual void drawPost();
	
public:
		
	MAKE_PARTICLESTUFF_NAME();
	
	// default constructor
	ShaderSilhouetteContour(Particles *ps=NULL);
	
	~ShaderSilhouetteContour(){	};
	
	void attachAttributes();

	
private:
		
	struct CPparams p;
};

#endif

