/*
@file ShaderMeshRenderer
@author Elmar Eisemann
This class draws an open mesh file.
It does not have to be triangulated!
*/

#ifndef PARTICLESHADERMESHRENDERER_H
#define PARTICLESHADERMESHRENDERER_H

#include "ParticleShader.h"
class MeshInterrogator;
class SurfaceMesh;

class ShaderMeshRenderer:public ParticleShader
{
public:
	MeshInterrogator * _mi;
	GLuint _glDisplayList;
	GLuint _glWireframeDisplayList;

	SurfaceMesh * _mesh;

	bool _wireframe;
	bool _lighting;

public:

	MAKE_PARTICLESTUFF_NAME();
	
	ShaderMeshRenderer(Particles *ps=NULL);

	virtual ~ShaderMeshRenderer();

	virtual void drawPost();
	void updateDisplayLists();
};

#endif
