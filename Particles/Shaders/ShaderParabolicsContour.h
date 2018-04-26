/**
 * Declaration of a shader for parabolic contours
 * @file ShaderParabolicsContour.h
 * @date 12/9/2005
 * @author Matei N. Stroila
 * @remarks
 */


#ifndef SHADERPARABOLICSCONTOUR_H
#define SHADERPARABOLICSCONTOUR_H


#include "ShaderContour.h"

class ShaderParabolicsContour : public ShaderContour
{	
protected:	
	virtual void findTangent(const unsigned int i);
	virtual void checkCriticalPoint4D(double *pstatus){};
	virtual void checkCriticalPoint5D(double *pstatus){};
	virtual void checkCriticalPointPC(){};
	virtual void drawPost();
	
public:
		
		MAKE_PARTICLESTUFF_NAME();
	
	// default constructor
	ShaderParabolicsContour(Particles *ps=NULL);
	
	~ShaderParabolicsContour();
	
	void attachAttributes();
	
	
private:
		
	bool setImplicit();	
	bool ImplicitIsSet;
	Implicit *imp1, *imp2;
};

#endif

