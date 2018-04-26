#include "ParticleShader.h"
#include "ParticleShaderChain.h"
#include "ImplicitInterrogator.h"

REGISTER_PARTICLESTUFF(ParticleShaderChain,"Shader:ParticleShaderChain");

ParticleShaderChain::ParticleShaderChain(Particles *ps)
	:ParticleShader(ps,std::string("ParticleShaderChain"))
{
	lineWidth=1;
	percentLength=0.1;
	autoUpdate=true;
}

int ParticleShaderChain::qlen()
{
	return 3;
}

void ParticleShaderChain::getq(double *q)
{
	q[0] = lineWidth;
	q[1] = percentLength;
	q[2] = autoUpdate;
}

void ParticleShaderChain::setq(double *q)
{
	lineWidth = q[0];
	percentLength = q[1];
	autoUpdate = (bool)q[2];
}

void ParticleShaderChain::qname(char **qn)
{
	qn[0] = "Line Width";
	qn[1] = "Limit for a chain to be removed";
	qn[2] = "Recalculated chain";
}

void ParticleShaderChain::draw(int s)
{
	// do not do this if too many particles exist
	if (ps->size()>500)
		autoUpdate=false;
	if (autoUpdate)
		findChains();
	glLineWidth(lineWidth);
	glDisable(GL_LIGHTING);
	// draw all lines
	glBegin(GL_LINES);
	for(unsigned int i=0;i<chains.size();i++)
	{
		gmVector3 vi=position->getPosition(i);
		gmVector3 v1=position->getPosition(chains[i].n1);
		gmVector3 v2=position->getPosition(chains[i].n2);
		glVertex3d(v1[0],v1[1],v1[2]);
		glVertex3d(vi[0],vi[1],vi[2]);
		glVertex3d(v2[0],v2[1],v2[2]);
		glVertex3d(vi[0],vi[1],vi[2]);
	}
	glEnd();
	glEnable(GL_LIGHTING);
}

void ParticleShaderChain::findChains()
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
