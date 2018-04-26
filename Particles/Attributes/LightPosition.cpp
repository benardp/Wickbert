/**
* Implementation of the light position attribute.
* @file LightPosition.cpp
* @date 27 May. 2005
* @author Matei N. Stroila
*/

#include "LightPosition.h"
#include "cleangl.h"

REGISTER_PARTICLESTUFF(LightPosition,"Attribute:LightPosition");

LightPosition::LightPosition(Particles *ps, const std::string& name):ParticleAttribute(ps,name)
{
	
	glGetLightfv(GL_LIGHT0, GL_POSITION, _lightPositionf);
	for(size_t i = 0; i < 4; i++) _lightPositionInput[i] = _lightPositionf[i];
	
	new PSParamgmVector4(this,&_lightPositionInput,_lightPositionInput,
		"light","light position","Light position in eye coordinate space.");
	new PSParamButton(this,new LightPositionSet(this),"set","Set Light Position", "Set Light Position");
	new PSParamBool(this,&_lightIsInWorldCoords,false,"worldCoord","isWorldCoords","The light will be placed in world coordinates, realize" 
																					"that the opengl light will remain at a relative position!"
																					"To update its position one needs to click on set light again");

} 

void LightPosition::prepare()
{
	if (!_lightIsInWorldCoords)
	{
		glGetLightfv(GL_LIGHT0, GL_POSITION, _lightPositionf);
		gmVector4 eyeLP;
		for(size_t i = 0; i < 4; i++) 
			eyeLP[i] = _lightPositionf[i];
		
		GLdouble modelview[16];
		glGetDoublev(GL_MODELVIEW_MATRIX,modelview);
		gmMatrix4 _m(modelview[0], modelview[4], modelview[8], modelview[12],
					modelview[1], modelview[5], modelview[9], modelview[13],
					modelview[2], modelview[6], modelview[10], modelview[14],
					modelview[3], modelview[7], modelview[11], modelview[15]);
		gmMatrix4 _m_inv = _m.inverse();
			
		_lightPosition = _m_inv * eyeLP;
	}
	else
	{
		_lightPosition = _lightPositionInput;
	}
}

void const LightPosition::getLightPosition(gmVector4& worldLightPosition)
{	
	worldLightPosition=_lightPosition;
} 

//set GL light to light position
//and thus also the actual light that we give back
void LightPosition::setLightPosition(void)
{
	if (!_lightIsInWorldCoords)
	{
		glMatrixMode(GL_MODELVIEW_MATRIX);
		glPushMatrix();
		glLoadIdentity();
		for(size_t i = 0; i < 4; i++) 
			_lightPositionf[i] = (float) _lightPositionInput[i];

		//we do not want to be ambiguous about directional and not directional lights
		if (gmIsZero(_lightPositionf[3]))
			_lightPositionf[3]=0;

		glLightfv(GL_LIGHT0, GL_POSITION, _lightPositionf);

		glPopMatrix();
	}
	else
	{
		for(size_t i = 0; i < 4; i++) 
			_lightPositionf[i] = (float) _lightPositionInput[i];

		//we do not want to be ambiguous about directional and not directional lights
		if (gmIsZero(_lightPositionf[3]))
			_lightPositionf[3]=0;

		glLightfv(GL_LIGHT0, GL_POSITION, _lightPositionf);		
	}
}

LightPosition::~LightPosition(void)
{
}
