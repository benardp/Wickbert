#include "ShowCell.h"
#include "ParticleLocalityGrid.h"

REGISTER_PARTICLESTUFF(ShowCell,"Shader:ShowCell");

ShowCell::ShowCell(Particles *ps)
	:ParticleShader(ps,std::string("ShowCell"))
{
	new Attached<ParticleLocalityGrid>(this,&plg);
}

void ShowCell::drawParticle(int i)
{
	if (!plg || i != ps->selectedParticle) return;

	gmVector3 a,b;

	plg->cellBox(i,a,b);

	glBegin(GL_LINES);
		glVertex3d(a[0],a[1],a[2]);
		glVertex3d(a[0],a[1],b[2]);

		glVertex3d(a[0],a[1],a[2]);
		glVertex3d(a[0],b[1],a[2]);

		glVertex3d(a[0],a[1],a[2]);
		glVertex3d(b[0],a[1],a[2]);

		glVertex3d(a[0],a[1],b[2]);
		glVertex3d(a[0],b[1],b[2]);

		glVertex3d(a[0],a[1],b[2]);
		glVertex3d(b[0],a[1],b[2]);

		glVertex3d(a[0],b[1],a[2]);
		glVertex3d(a[0],b[1],b[2]);

		glVertex3d(a[0],b[1],a[2]);
		glVertex3d(b[0],b[1],a[2]);

		glVertex3d(b[0],a[1],a[2]);
		glVertex3d(b[0],a[1],b[2]);

		glVertex3d(b[0],a[1],a[2]);
		glVertex3d(b[0],b[1],a[2]);

		glVertex3d(a[0],b[1],b[2]);
		glVertex3d(b[0],b[1],b[2]);

		glVertex3d(b[0],a[1],b[2]);
		glVertex3d(b[0],b[1],b[2]);
		
		glVertex3d(b[0],b[1],a[2]);
		glVertex3d(b[0],b[1],b[2]);
	glEnd();
}
