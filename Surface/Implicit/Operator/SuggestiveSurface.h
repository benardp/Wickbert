/*
 *  SuggestiveSurface.h
 *  gps
 *
 *  Created by Matei Stroila on 6/20/05.
 *  
 *
 */

#ifndef SUGGESTIVESURFACE_H
#define SUGGESTIVESURFACE_H

class SuggestiveSurface : public Implicit
{
public:
	SuggestiveSurface();
	SuggestiveSurface(Implicit *F, gmVector3 *cameraP); ///< Explicit constructor.
	
	double proc(const gmVector3 & x);
	Intervald proc(const Box<double>& x);
	gmVector3 grad(const gmVector3 & x);
	
	~SuggestiveSurface(){}; 
	
	Particles* p;
	Implicit *myF; ///< The input implicit.
	
	//	void getq(double*);
	void _setq(double*);
	//	int qlen() { return 3; };
	//	void getqname(char** qn);
	
	MAKE_NAME();
	
	
	
private:
		
		gmVector3 *myV; ///< The input view direction.
	ViewDependence *view;
	
	
};

#endif

