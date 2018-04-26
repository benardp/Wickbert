#include "IntervalSubdivision.h"
#include "ImplicitInterrogator.h"
#include "ParticleBoundingBox.h"

REGISTER_PARTICLESTUFF(IntervalSubdivision,"Shader:IntervalSubdivision");

IntervalSubdivision::IntervalSubdivision(Particles *ps)
	:ParticleShader(ps,std::string("IntervalSubdivision"))
{
	new PSParamDouble(this,&res,1.0,"res","subdiv. res.",
		"Resolution below which intervals are no longer subdivided.");
	new PSParamgmVector3(this,&color,gmVector3(1.0,0.0,0.0),"color","box color",
		"Color of the displayed box.");
	new Attached<ParticleBoundingBox>(this,&pbox);
	new Attached<ImplicitInterrogator>(this,&imp_int_A,"Implicit A","impA","Implicit A",
		"One of the implicits on which to detect a zero.");
	new Attached<ImplicitInterrogator>(this,&imp_int_B,"Implicit B","impB","Implicit B",
		"One of the implicits on which to detect a zero.");
	new Attached<ImplicitInterrogator>(this,&imp_int_C,"Implicit C","impC","Implicit C",
		"One of the implicits on which to detect a zero.");
	new Attached<ImplicitInterrogator>(this,&imp_int_D,"Implicit D","impD","Implicit D",
		"One of the implicits on which to detect a zero.");
}

void IntervalSubdivision::drawPre()
{
	impA = imp_int_A ? imp_int_A->getImplicit() : NULL;
	impB = imp_int_B ? imp_int_B->getImplicit() : NULL;
	impC = imp_int_C ? imp_int_C->getImplicit() : NULL;
	impD = imp_int_D ? imp_int_D->getImplicit() : NULL;

	if (impB == impA) impB = NULL;
	if (impC == impB || impC == impA) impC = NULL;
	if (impD == impC || impD == impB || impD == impA) impD = NULL;

	Box3d b;
	if (pbox) b = Box3d(pbox->min,pbox->max);

	glDisable(GL_LIGHTING);
	glColor3d(color[0],color[1],color[2]);

	subdivide(b);

	glEnable(GL_LIGHTING);
}

void IntervalSubdivision::subdivide(Box3d b)
{
	gmVector3 l = b.low();
	gmVector3 h = b.high();

	// draw wireframe box

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

//	if (b.width() > res &&

	double b0 = b[0].width();
	double b1 = b[1].width();
	double b2 = b[2].width();

	if ((b0 > res || b1 > res || b2 > res) &&
		(impA == NULL || impA->proc(b).contains(0.0)) &&
		(impB == NULL || impB->proc(b).contains(0.0)) &&
		(impC == NULL || impC->proc(b).contains(0.0)) &&
		(impD == NULL || impD->proc(b).contains(0.0))) {

		std::vector<Box3d> boxes;
		int i,n;
		gmVector3 c = 0.5*(l + h);

		if (b0 > res) {
			boxes.push_back(Box3d(gmVector3(l[0],l[1],l[2]),gmVector3(c[0],h[1],h[2])));
			boxes.push_back(Box3d(gmVector3(c[0],l[1],l[2]),gmVector3(h[0],h[1],h[2])));
		} else {
			boxes.push_back(b);
		}
		if (b1 > res) {
			n = boxes.size();
			for (i = 0; i < n; i++) {
				Box3d a = boxes[i];
				boxes[i] = Box3d(a.low(),gmVector3(a[0].high(),c[1],a[2].high()));
				boxes.push_back(Box3d(gmVector3(a[0].low(),c[1],a[2].low()),a.high()));
			}
		}
		if (b2 > res) {
			n = boxes.size();
			for (i = 0; i < n; i++) {
				Box3d a = boxes[i];
				boxes[i] = Box3d(a.low(),gmVector3(a[0].high(),a[1].high(),c[2]));
				boxes.push_back(Box3d(gmVector3(a[0].low(),a[1].low(),c[2]),a.high()));
			}
		}

		n = boxes.size();
		for (i = 0; i < n; i++) {
			subdivide(boxes[i]);
		}
	}
}