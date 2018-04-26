/**
 * Implementation of a new shader that draws the bounding box of a Particles object
 * @file ShaderBoundingBox.cpp
 * @date 03/22/2006
 * @author Matei Stroila
 * @remarks
 */

#include "ShaderBoundingBox.h"
#include "ParticleBoundingBox.h"
#include "Surface/Box.h"

REGISTER_PARTICLESTUFF(ShaderBoundingBox,"Shader:ShaderBoundingBox");

ShaderBoundingBox::ShaderBoundingBox(Particles *ps)
	:ParticleShader(ps,std::string("ShaderBoundingBox"))
{
	new PSParamgmVector3(this,&color,gmVector3(1.0,0.0,0.0),"color","box color",
		"Color of the displayed box.");
	new Attached<ParticleBoundingBox>(this,&pbox);
}


void ShaderBoundingBox::drawPre()
{
	if (!pbox) return; 
	Box3d b(pbox->min,pbox->max);
	gmVector3 l(b.low());
	gmVector3 h(b.high());

	// draw wireframe box
	glDisable(GL_LIGHTING);
	glColor3d(color[0],color[1],color[2]);
	glBegin(GL_LINE_LOOP);
		glVertex3d(l[0],l[1],l[2]);
		glVertex3d(h[0],l[1],l[2]);
		glVertex3d(h[0],h[1],l[2]);
		glVertex3d(l[0],h[1],l[2]);
	glEnd();

	glBegin(GL_LINES);
		glVertex3d(l[0],l[1],l[2]);	glVertex3d(l[0],l[1],h[2]);
		glVertex3d(h[0],l[1],l[2]);	glVertex3d(h[0],l[1],h[2]);
		glVertex3d(h[0],h[1],l[2]);	glVertex3d(h[0],h[1],h[2]);
		glVertex3d(l[0],h[1],l[2]);	glVertex3d(l[0],h[1],h[2]);
	glEnd();

	glBegin(GL_LINE_LOOP);
		glVertex3d(l[0],l[1],h[2]);
		glVertex3d(h[0],l[1],h[2]);
		glVertex3d(h[0],h[1],h[2]);
		glVertex3d(l[0],h[1],h[2]);
	glEnd();
	glEnable(GL_LIGHTING);

}
