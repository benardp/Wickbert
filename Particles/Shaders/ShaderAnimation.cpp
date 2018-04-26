/**
* Implementation of the canvas animation.
* @file ShaderAnimation.cpp
* @date 14 November. 2005
* @author Matei N. Stroila
*/


#include "ShaderAnimation.h"

REGISTER_PARTICLESTUFF(ShaderAnimation,"Shader:ShaderAnimation");

ShaderAnimation::ShaderAnimation(Particles *ps)
:ParticleShader(ps,std::string("ShaderAnimation"))
{
	new PSParamBool(this,&rotateX,true,"rotateX","Rotate about the x-axis","Set true if want to rotate the camera about the x-axis");
	new PSParamBool(this,&rotateY,false,"rotateY","Rotate about the y-axis","Set true if want to rotate the camera about the y-axis");
	new PSParamBool(this,&rotateZ,false,"rotateZ","Rotate about the z-axis","Set true if want to rotate the camera about the z-axis");
	new PSParamDouble(this,&speed,1.0,"speed","Rotation speed","Set the rotation speed");
	new PSParamDouble(this,&zoom,-5,"zoom","Camera zoom","Set the zoom");
	new PSParamButton(this,new ShaderAnimationStart(this),"start","Start Animation", "Start the animation.");
	new PSParamButton(this,new ShaderAnimationStop(this),"stop","Stop Animation",
					  "Stop the animation.");
}


void ShaderAnimation::startAnimation(){
	ps->particleSystem->getEulerAnglesAndZoom(&xinit,&yinit,&zoom);
	ps->particleSystem->particleSystems->setAnimation(true, rotateX, rotateY, rotateZ, speed, zoom, xinit, yinit);	
}

void ShaderAnimation::stopAnimation(){
	ps->particleSystem->getEulerAnglesAndZoom(&xinit,&yinit,&zoom);
	ps->particleSystem->particleSystems->setAnimation(false, rotateX, rotateY, rotateZ, speed, zoom, xinit, yinit);	
}
