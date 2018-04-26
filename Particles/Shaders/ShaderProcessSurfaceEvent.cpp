/*
 *  ShaderProcessSurfaceEvent.cpp
 *  GPS
 *
 *  Created by Matei Stroila on 11/16/05.
 *  
 *
 */

#include "ShaderProcessSurfaceEvent.h"

REGISTER_PARTICLESTUFF(ShaderProcessSurfaceEvent,"Shader:ShaderProcessSurfaceEvent");

ShaderProcessSurfaceEvent::ShaderProcessSurfaceEvent(Particles *ps)
:ParticleShader(ps,std::string("ShaderProcessSurfaceEvent"))
{
	new PSParamBool(this,&processFlag,false,"process","Process Event","Set true if want to process some event in Surface");
	new PSParamString(this,&stringEvent,"<empty>","stringEvent","String passed", "String passed to Surfaces.");
	sfcs = NULL;
}


void ShaderProcessSurfaceEvent::drawPost()
{
	if(!processFlag) return;
	sfcs = ps->particleSystem->surfaces;
	Surfaces::iterator iter;
	for(iter = (*sfcs).begin(); iter != (*sfcs).end(); ++iter)
	{
		(*iter)->processEvent(stringEvent);
	}
	processFlag = false;
}