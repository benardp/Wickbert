#include "ParticleShader.h"
#include "ParticleShaderStroke.h"
#include "SilhouetteAdhesion.h"

REGISTER_PARTICLESTUFF(ParticleShaderStroke,"Shader:ParticleShaderStroke");

ParticleShaderStroke::ParticleShaderStroke(Particles *ps)
	:ParticleShader(ps,std::string("ParticleShaderStroke"))
{
	lineWidth=0.01;
	percentLength=0.1;
	autoUpdate=true;
	inkDarkness=0.5;
	curve=0.001;
	lineTaper=0.01;
}

int ParticleShaderStroke::qlen()
{
	return 6;
}

void ParticleShaderStroke::getq(double *q)
{
	q[0] = lineWidth;
	q[1] = percentLength;
	q[2] = autoUpdate;
	q[3] = curve;
	q[4] = inkDarkness;
	q[5] = lineTaper;
}

void ParticleShaderStroke::setq(double *q)
{
	lineWidth = q[0];
	percentLength = q[1];
	autoUpdate = (bool)q[2];
	curve = q[3];
	inkDarkness = q[4];
	lineTaper = q[5];
}

void ParticleShaderStroke::qname(char **qn)
{
	qn[0] = "Line Width";
	qn[1] = "Limit for a chain to be removed";
	qn[2] = "Recalculated chain";
	qn[3] = "Curve";
	qn[4] = "Ink Darkness";
	qn[5] = "Line Taper";
}

void ParticleShaderStroke::draw(int s)
{
	// do not do this if too many particles exist
	if (ps->size()>500)
		autoUpdate=false;
	if (autoUpdate)
		findChains();


	glDisable(GL_LIGHTING);
	// draw all lines
	for(unsigned int i=0;i<chains.size();i++)
	{
		gmVector3 ni=orient->getNormal(i);
		gmVector3 vi=position->getPosition(i);
		gmVector3 v1=position->getPosition(chains[i].n1);
		gmVector3 v2=position->getPosition(chains[i].n2);
		gmVector3 v3=(vi+v1+v1)*0.3333333333333333333;
		gmVector3 v4=(vi+v2+v2)*0.3333333333333333333;
		gmVector3 v5=(vi+v3)*0.5+ni*v3.length()*curve;
		gmVector3 v6=(vi+v4)*0.5+ni*v4.length()*curve;
		gmVector3 lineDir=(v1-vi)*0.2;
		// outside line
		gmVector3 vv3=v3+ni*lineWidth+ni*lineTaper*2;
		gmVector3 vv5=v5+ni*lineWidth+ni*lineTaper*1.5;
		gmVector3 vvi=vi+ni*lineWidth+ni*lineTaper;
		gmVector3 vv6=v6+ni*lineWidth+ni*lineTaper*0.5;
		gmVector3 vv4=v4+ni*lineWidth+ni*lineTaper*0.2;

		float color[3];
		color[0]=0.0;
		color[1]=0.0;
		color[2]=inkDarkness;

		glBegin(GL_QUADS);
		glColor3fv(color);
		glVertex3f(vv3[0],vv3[1],vv3[2]);
		glVertex3f(v3[0],v3[1],v3[2]);
		color[2]=color[2]+0.1;
		glColor3fv(color);
		glVertex3f(v5[0],v5[1],v5[2]);
		glVertex3f(vv5[0],vv5[1],vv5[2]);

		glVertex3f(vv5[0],vv5[1],vv5[2]);
		glVertex3f(v5[0],v5[1],v5[2]);
		color[2]=color[2]+0.1;
		glColor3fv(color);
		glVertex3f(vi[0],vi[1],vi[2]);
		glVertex3f(vvi[0],vvi[1],vvi[2]);

		glVertex3f(vvi[0],vvi[1],vvi[2]);
		glVertex3f(vi[0],vi[1],vi[2]);
		color[2]=color[2]+0.1;
		glColor3fv(color);
		glVertex3f(v6[0],v6[1],v6[2]);
		glVertex3f(vv6[0],vv6[1],vv6[2]);

		glVertex3f(vv6[0],vv6[1],vv6[2]);
		glVertex3f(v6[0],v6[1],v6[2]);
		color[2]=color[2]+0.1;
		glColor3fv(color);
		glVertex3f(v4[0],v4[1],v4[2]);
		glVertex3f(vv4[0],vv4[1],vv4[2]);
		glEnd();
	}
	glEnable(GL_LIGHTING);
}

void ParticleShaderStroke::findChains()
{
	unsigned int i,j;
	double averageLength=0;

	chains.clear();
	// find the closest two
	for(i=0;i<ps->size();i++)
	{
		Neighbor n;
		n.n1=0;
		n.n2=0;
		// set minDistance to max
		double minDistance;
		double curDistance;
		// smallest
		minDistance=gmGOOGOL;
		for(j=0;j<ps->size();j++)
		{
			// not the same particle
			if (i!=j)
			{
				curDistance=distance(position->getPosition(j),position->getPosition(i));
				if (curDistance<minDistance)
				{
					n.n1=j;
					minDistance=curDistance;
				}
			}
		}
		averageLength+=minDistance;

		// second smallest
		minDistance=gmGOOGOL;
		for(j=0;j<ps->size();j++)
		{
			// not the same particles
			if (i!=j && j!=n.n1)
			{
				curDistance=distance(position->getPosition(j),position->getPosition(i));
				if (curDistance<minDistance)
				{
					n.n2=j;
					minDistance=curDistance;
				}
			}
		}
		averageLength+=minDistance;
		chains.push_back(n);
	}

	// find average
	averageLength=averageLength/chains.size();
	double limitLength=averageLength/2.0;
	for(i=0;i<chains.size();i++)
	{
		// find long links and remove them by link to self
		double d;
		d=distance(position->getPosition(i),position->getPosition(chains[i].n1))+distance(position->getPosition(i),position->getPosition(chains[i].n1));
		d=d*percentLength;
		if (d>limitLength)
		{
			chains[i].n1=i;
			chains[i].n2=i;
		}
	}
}
