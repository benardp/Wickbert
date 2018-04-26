/*
 *  RotatedImplicit.h
 *  gps
 *
 *  Created by Matei Stroila on 10/13/05.
 *  
 */

#ifndef RotatedImplicit_h
#define RotatedImplicit_h

#include "Surface/Implicit/Implicit.h"

#if defined(_MSC_VER)  // Identifying MS VC++ compilers
         // Special headers for MS compilers  
	#include <wtypes.h>
	#include <wingdi.h>  

#endif
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif


class RotatedImplicit : public Implicit
{
public:
	RotatedImplicit();
	RotatedImplicit(Implicit *F, gmMatrix3& R); ///< Explicit constructor.
	
	virtual double proc(const gmVector3 & x);
	virtual Intervald proc(const Box<double>&  b);
	virtual gmVector3 grad(const gmVector3 & x);
	virtual gmMatrix3 hess(const gmVector3 & x);
	
	void setF(Implicit *F);
	void setRotation(gmMatrix3& R);
	
	~RotatedImplicit(){}; 
	MAKE_NAME();

 private:
		
	Implicit *myF; ///< The input implicit.
	gmMatrix3 myR; ///< Rotation Matrix
	bool rotHasChanged;
	void checkRotation(void);
};

#endif

