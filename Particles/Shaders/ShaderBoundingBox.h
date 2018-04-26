/**
 * Declaration of a new shader that draws the bounding box of a Particles object
 * @file ShaderBoundingBox.h
 * @date 03/22/2006
 * @author Matei Stroila
 * @remarks
 */

#ifndef ShaderBoundingBox_h
#define ShaderBoundingBox_h

#include "ParticleShader.h"
class ParticleBoundingBox;

class ShaderBoundingBox : public ParticleShader
{
public:
	MAKE_PARTICLESTUFF_NAME();

	ShaderBoundingBox(Particles *ps=NULL);

	virtual void drawPre();

private:

	ParticleBoundingBox *pbox;
	gmVector3 color;

};

#endif
