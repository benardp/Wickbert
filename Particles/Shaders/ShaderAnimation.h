/**
* Declaration of the canvas animation.
* @file ShaderAnimation.h
* @date 14 November. 2005
* @author Matei N. Stroila
*/

#ifndef SHADERANIMATION_H
#define SHADERANIMATION_H

#include "ParticleShader.h"

class ShaderAnimation : public ParticleShader
{
	
public:
		
	MAKE_PARTICLESTUFF_NAME();
	
	ShaderAnimation(Particles *ps=NULL);
	~ShaderAnimation(){};
	
	void startAnimation();
	void stopAnimation();
	
private:
	
	bool rotateX;
	bool rotateY;
	bool rotateZ;	
	double speed;
	double zoom;
	double xinit,yinit;
	
};

class ShaderAnimationStart : public PSParamButton::Callback
{
public:
	ShaderAnimation *sha;
	ShaderAnimationStart(ShaderAnimation *me) {sha = me;}
	virtual void onbuttonpress() {sha->startAnimation();}
};

class ShaderAnimationStop : public PSParamButton::Callback
{
public:
	ShaderAnimation *sha;
	ShaderAnimationStop(ShaderAnimation *me) {sha = me;}
	virtual void onbuttonpress() {sha->stopAnimation();}
};
#endif