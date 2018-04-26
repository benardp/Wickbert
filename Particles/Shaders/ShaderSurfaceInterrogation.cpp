#include "ShaderSurfaceInterrogation.h"
#include "ImplicitInterrogator.h"
#include "ParticlePosition.h"

REGISTER_PARTICLESTUFF(ShaderSurfaceInterrogation,"Shader:ShaderSurfaceInterrogation");

ShaderSurfaceInterrogation::ShaderSurfaceInterrogation(Particles *ps)
	:ParticleShader(ps,std::string("ShaderSurfaceInterrogation"))
{
	new Attached<ImplicitInterrogator>(this,&_implicitInter);
	new Attached<ParticlePosition>(this,&_particlePosition);
	

	new PSParamBool(this, &_updateValuesAtEachIter,true,"update","Update particles", "Update the particles displayed at each iteration");
	new PSParamBool(this, &_activateBlend,true,"Blend","Blending", "Activate GL Blending");
	new PSParamBool(this, &_activateGradient,false,"grad","use grad", "Use line segments to visualize the gradient");
	new PSParamDouble(this, &_gradientScale,0.1,"gradScale","scale gradient", "Scale the line segments visualizing the gradient");
	new PSParamDouble(this, &_alphaScaling,1.0,"alpha","alpha scaling", "The implicit functions value (eventually after normalization) is scaled with this value");

	new PSParamDouble(this, &_size,1.0,"size","size", "PixelSize of the primitives");

	new PSParamDouble(this, &_upperThreshold,FLT_MAX,"hThres","high thresh.", "Upper threshold for culling of particles");
	new PSParamDouble(this, &_lowerThreshold,-FLT_MAX,"lThres","lower thresh.", "Lower threshold for culling of particles");
	new PSParamBool(this, &_normalizeValues,false,"normalize","normalize values", "Normalize the values between upper and lower threshold");


	new PSParamComboBox(this,new ValueSelectorCallback(this),"value","func value", "Select the value that is displayed by the particles");


	_glDisplayListIsUpToDate=false;
}

ShaderSurfaceInterrogation::~ShaderSurfaceInterrogation()
{
}

void ShaderSurfaceInterrogation::drawPre()
{
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glLineWidth((float)_size);
	glPointSize((float)_size);

	glDisable(GL_LIGHTING);

	if (_activateBlend)
	{
		glBlendFunc(GL_ZERO, GL_SRC_ALPHA);
		glEnable(GL_BLEND);
		glDisable(GL_DEPTH_TEST);
	}
	
	_tempThresDiff=_upperThreshold-_lowerThreshold;

	if (_updateValuesAtEachIter)
	{
		if (_glDisplayListIsUpToDate)
		{
			_glDisplayListIsUpToDate=false;	
		}
	}
	else
	{
		if (!_glDisplayListIsUpToDate)
		{
			//we create  a new display list, the old is not up to date
			glDeleteLists(1,_glDisplayList);
			_glDisplayList= glGenLists(1);
			glNewList(_glDisplayList,GL_COMPILE);
			//if we do not set this, the particles won't be drawn
			_updateValuesAtEachIter=true;
			for (unsigned int i=0;i<ps->size();++i)
				drawParticle(i);
			_updateValuesAtEachIter=false;
			glEndList();
			_glDisplayListIsUpToDate=true;
		}

		//we draw all the particles
		glCallList(_glDisplayList);
	}


}

void ShaderSurfaceInterrogation::drawPost()
{
	glPopAttrib();
}

void ShaderSurfaceInterrogation::drawParticle(int i)
{
	if (_particlePosition==0)
		return;

	//if we have a valid displaylist, we can safely bail out
	//this is why I would love to put this drawParticle function
	//also in the predraw... but anyway... let's just stick to the 
	//original program structure
	if (!_updateValuesAtEachIter)
		return;

	double proc=_implicitInter->proc(i);
	if (proc <_lowerThreshold)
		return;

	if (proc >_upperThreshold)
		return;
	
	if (_normalizeValues)
		proc=(proc-_lowerThreshold)/(_tempThresDiff);

	if (_activateBlend)
	{
		glColor4d(1,0,0,proc*_alphaScaling);

	}
	else
	{
		glColor3d(proc*_alphaScaling,0,0);
	}

	const gmVector3 pos=_particlePosition->getPosition(i);
	if (_activateGradient)
	{
		glBegin(GL_LINES);
			glVertex3d(pos[0],pos[1],pos[2]);
			gmVector3 gradient=_implicitInter->grad(i);
			glVertex3d(  pos[0]+_gradientScale*gradient[0]
						,pos[1]+_gradientScale*gradient[1]
						,pos[2]+_gradientScale*gradient[2]);
		glEnd();
	}
	else
	{
		glBegin(GL_POINTS);
			glVertex3d(pos[0],pos[1],pos[2]);
		glEnd();
	}

}



ShaderSurfaceInterrogation::ValueSelectorCallback::ValueSelectorCallback(ShaderSurfaceInterrogation * shader)
: PSParamComboBox::Callback()
, _shader(shader)
{
	//We could replace this by an own factory... but to me it seems overkill...
	_choices.push_back(std::string("proc"));
}

ShaderSurfaceInterrogation::ValueSelectorCallback::~ValueSelectorCallback()
{
}

void ShaderSurfaceInterrogation::ValueSelectorCallback::itemselected(unsigned int selection)
{
	// for the moment only proc is used.
}
