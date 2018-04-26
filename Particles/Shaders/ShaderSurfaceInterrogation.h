/*
@file ShaderITKVolumeData
@author Elmar Eisemann
This class draws an itk implicit as a point cloud using simple 
density based blending.
*/

#ifndef PARTICLESHADERSURFACEINTERROGATION_H
#define PARTICLESHADERSURFACEINTERROGATION_H

#include "ParticleShader.h"
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif
class ImplicitInterrogator;
class ParticlePosition;


class ShaderSurfaceInterrogation:public ParticleShader
{
protected:
	ImplicitInterrogator * _implicitInter;
	ParticlePosition * _particlePosition;

	//blend particle contributions together using OpenGL
	bool _activateBlend;
	double _alphaScaling;
	//represent gradient with line segments
	bool _activateGradient;
	double _gradientScale;

	double _size;

	//specify bounds for the particles to draw
	double _upperThreshold;
	double _lowerThreshold;
	bool _normalizeValues;

	//at each iteration the values of the particles change
	//otherwise just keep them
	bool _updateValuesAtEachIter;


	//temporary variables.
	double _tempThresDiff;

	GLuint _glDisplayList;
	bool _glDisplayListIsUpToDate;

class ValueSelectorCallback : public PSParamComboBox::Callback
{
protected:
	ShaderSurfaceInterrogation * _shader;
public:
	ValueSelectorCallback(ShaderSurfaceInterrogation * surface);
	virtual ~ValueSelectorCallback();
	virtual void itemselected(unsigned int i);
};

public:

	MAKE_PARTICLESTUFF_NAME();
	
	ShaderSurfaceInterrogation(Particles *ps=NULL);
	virtual ~ShaderSurfaceInterrogation();

	virtual void drawParticle(int i);
	virtual void drawPre();
	virtual void drawPost();
};
#endif
