/*
 *  ShaderProcessSurfaceEvent.h
 *  GPS
 *
 *  Created by Matei Stroila on 11/16/05.
 *  
 *
 */

#ifndef SHADERPROCESSSURFACEEVENT_H
#define SHADERPROCESSSURFACEEVENT_H

#include "ParticleShader.h"
#include "ParticleSystem.h"
#include "Surface/Surface.h"

class ShaderProcessSurfaceEvent : public ParticleShader
{
	
public:
	
	MAKE_PARTICLESTUFF_NAME();
	
	ShaderProcessSurfaceEvent(Particles *ps=NULL);
	~ShaderProcessSurfaceEvent(){};
	
	virtual void drawPost();
	
private:
		
	bool processFlag; //flag for event capture
	Surfaces *sfcs;
	std::string stringEvent;
	
	
};

#endif