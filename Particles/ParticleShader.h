/**
 * @file ParticleShader.h
 * @author Wen Su
 */
#ifndef PARTICLESHADER_H
#define PARTICLESHADER_H

#ifdef USE_CLEANGL
#include "cleangl.h"
#else
#if defined(_WIN32)
#include <windows.h>
#endif
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif
#endif

#include "Particles.h"
#include "ParticleStuff.h"

enum {
	WB_CLICK = 1,
	WB_DOUBLECLICK,
	WB_DRAG,
	WB_RIGHTCLICK,
	WB_SHIFT = 8,
	WB_CTRL = 16,
	WB_ALT = 32,
};


/** ParticleShader describes how particles interface with the user,
 * either through visual output or through event driven input.
 *
 * The Particles object draws particles using the following control flow:\n
 \verbatim
 foreach shader s
   s->drawPre()

 foreach particle i
   foreach shader s
     s->drawParticle(i)

 foreach shader s
   s->drawPost()
 \endverbatim
 * 
 * User input is controlled by specializations of ParticleShader that redefine
 * the event method.
 */
class ParticleShader : public ParticleStuff
{

public:

	/// default constructor
	ParticleShader(Particles *ps=NULL,const std::string& name=std::string("ParticleShader"));

	virtual void attachAttributes(); 

	virtual ~ParticleShader() {}

	/// Add shader to a particle system.
	virtual void addMeTo(Particles *new_ps);

	virtual void setParticleSystem(Particles *new_ps);

	/// Remove the shader from the current particle system.
	virtual void removeMe();

	virtual bool moveAfter(ParticleStuff *stuff);

	/** Process an event. Defaults to a no-op.
	 */
	virtual void event(int e) { }

	/** Called in shader order before any drawParticle's
	 */
	virtual void drawPre() {}

	/** Called in particle order in shader order
	 */
	virtual void drawParticle(int i) {}

	/** Called in shader order after all drawParticle's are finished
	 */
	virtual void drawPost() {}

	virtual void clear() {}

	int iterations;
	double pre_time, draw_time, post_time;
};


#endif
