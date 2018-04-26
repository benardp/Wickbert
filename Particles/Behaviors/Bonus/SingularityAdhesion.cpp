#include <cmath>
// this is a bad implementation it depends on many shaders, it should not
#include "ParticlesInclude.h"

REGISTER_PARTICLESTUFF(SingularityAdhesion,"Behavior:SingularityAdhesion");

// default constructor
SingularityAdhesion::SingularityAdhesion(Particles *ps)
: ParticleBehavior(ps, std::string("SingularityAdhesion"))
{
	// set default param
	movTolerance=0.1;
	popDifTolerance=3;
	maxSingularPopulation=20;
	stableTolerance=30;
	recalculate=0;
	distnaceTol=0.1;
	normalTol=0.1;
	isStable=false;
	impInt=NULL;
	target=NULL;
}

void SingularityAdhesion::attachAttributes()
{
	ParticleBehavior::attachAttributes();
	attachAttribute(p_orient,std::string("ParticleOrientation"));
	attachAttribute(impInt,std::string("ImplicitInterrogator"));
}

void SingularityAdhesion::createSingularPS(Particles* ps)
{
	// there is a target or currect ps is null return
	if (ps==NULL || target!=NULL)
		return;
	// let the floater knows the control
	target=new Particles(std::string("Singularity"));
	new ParticleShaderSphere(target);
    // find particle repulsion
	ParticleRepulsion *particleRpl=NULL;
	particleRpl=ps->getBehavior<ParticleRepulsion>(std::string("ParticleRepulsion"));
	SingularityRepulsion *sinRel=new SingularityRepulsion(target, ps, particleRpl);
	sinRel->setParticleSystem(target);
	ImplicitInterrogator *tempImpInt=new ImplicitInterrogator(target, impInt->getImplicit());
	tempImpInt->setParticleSystem(target);
	target->createParticleStuff<SingularityInterrogator>();
	target->createParticleStuff<SurfaceAdhesion>();
	target->createParticleStuff<VelocityIntegrator>();
	//target->createParticleStuff<ParticleFate>();

	// use this to find neighbors.
	sinChn=target->createParticleStuff<ParticleShaderChain>();
	target->createParticleStuff<ParticleShaderSphere>();

	// record it in system
	if (ps->particleSystem)
		ps->particleSystem->particles.push_back(target);
}

/// parameters
int SingularityAdhesion::qlen()
{
	return 7;
}

void SingularityAdhesion::getq(double *q)
{
	q[0]=distnaceTol;
	q[1]=normalTol;
	q[2]=popDifTolerance;
	q[3]=movTolerance;
	q[4]=maxSingularPopulation;
	q[5]=stableTolerance;
	q[6]=recalculate;
}

void SingularityAdhesion::setq(double *q)
{
	distnaceTol=q[0];
	normalTol=q[1];
	popDifTolerance=q[2];
	movTolerance=q[3];
	maxSingularPopulation=q[4];
	stableTolerance=q[5];
	recalculate=q[6];
}

void SingularityAdhesion::qname(char **qn)
{
	qn[0]="distance traveled tolerance";
	qn[1]="angels in radient tolerance";
	qn[2]="stable population difference tolerance";
	qn[3]="movement tolerance";
	qn[4]="max singularities population";
	qn[5]="stable tolerance";
	qn[6]="set 1 to recalculate";
}

void SingularityAdhesion::applyForce()
{
	// make sure the history is the same size as the particle system
	updateHistorySize();

	// reset to recalculate singular points
	if (recalculate)
	{
		createSingularPS(ps);
		resetSingularities();
	}
	// do not do anything unless it is stable
	if (!isStable)
		return;

	oldN.resize(ps->size());
	// for each particle just record the position
	for (unsigned int i = 0; i < ps->size(); i++)
		// record the normal before it is updated
		oldN[i]=p_orient->getNormal(i);
}

void SingularityAdhesion::applyConstraint()
{
	// for each particle just record the position
	double q[2];
	movement=0;
	for (unsigned int i = 0; i < ps->size(); i++)
	{
		getqpp(q,i);
		// add the distance traveled to half of the old value
		double templ=velocity->v[i].length();
		movement+=templ;
		q[0]=velocity->v[i].length()+q[0]*0.5;
		// set the param
		setqpp(q,i);
	}
}

void SingularityAdhesion::cleanup()
{
	// do not do anything unless it is stable
	if (!isStable)
		return;

	// remove old points
	if (target==NULL)
		return;

	double q[2];
	double angle;
	unsigned int i;

	int oldSize=singularIndex.size();
	singularPoints.resize(ps->size());
	singularIndex.resize(ps->size());
	for (i = oldSize; i < (int)singularIndex.size(); i++)
	{
		singularIndex[i]=0;
	}

	for (i = 0; i < ps->size(); i++)
	{
		// for each particle just record the position
		getqpp(q,i);
		// record the new normal
		angle=acos(dot(oldN[i],p_orient->getNormal(i)));
		// add the angle moved
		q[1]=(float)angle+q[1]*0.5;
		// set the param
		setqpp(q,i);
		
		ParticlePosition *targetPosition=target->getAttribute<ParticlePosition>(std::string("ParticlePosition"));
		// check if it is in stable state and not exceed max population
		if (isStable && (int)target->size()<maxSingularPopulation)
		{
			// add singular points to the control particle system.
			if (distance[i] > distnaceTol || normalChange[i] > normalTol)
			{
				singularIndex[i]=singularIndex[i]+1;
				int temp=singularIndex[i];
				if (singularIndex[i]>stableTolerance)
				{
					// first timer use the value
					if (singularIndex[i]==(stableTolerance+1))
					{
						singularPoints[i]=position->getPosition(i);
						target->addParticle();
						targetPosition->setPosition(target->size()-1,singularPoints[i]);
						// reset to create new singular points
						singularIndex[i]=0;

						// max population reached
						if (target->size()==maxSingularPopulation)
							sinChn->findChains();
					}
					// else use the average
					else
					{
						singularPoints[i]=((singularIndex[i]-1)*singularPoints[i]+position->getPosition(i))/(float)singularIndex[i];
					}
					//// add this particle only if this index has exceeded the tolerance for the first time
					//if (singularIndex[i]==stableTolerance+1)
					// clamp the index
				}
			}
		}
	}
}

void SingularityAdhesion::resetSingularities()
{
	// reset stable flag
	isStable=true;
	// reset recalculate flag
	recalculate=0;
	// clear points
	if (target)
		target->clear();
	// reset history
	singularPoints.clear();
	singularIndex.clear();
}

/// per particle data
int SingularityAdhesion::qlenpp()
{
	return 2;
}

void SingularityAdhesion::getqpp(double *q, int i)
{
	q[0] = distance[i];
	q[1] = normalChange[i];
}
void SingularityAdhesion::setqpp(double *q, int i)
{
	distance[i] = q[0];
	normalChange[i] = q[1];
}

void SingularityAdhesion::qnamepp(char **qn)
{
	qn[0] = "distance traveled over last n iteration";
	qn[1] = "normal angle change in radient";
}

void SingularityAdhesion::particleAdded()
{
	distance.push_back(0);
	normalChange.push_back(0);
}
	
void SingularityAdhesion::particleRemoved(unsigned int i)
{
	distance[i] = distance.back();
	distance.pop_back();

	normalChange[i] = normalChange.back();
	normalChange.pop_back();
}

void SingularityAdhesion::updateHistorySize()
{
	unsigned int oldSize=distance.size();
	unsigned int newSize=ps->size();
	if (oldSize!=newSize)
	{
		distance.resize(newSize);
		normalChange.resize(newSize);
		for(unsigned int i=oldSize;i<newSize;i++)
		{
			distance[i]=0.0;
			normalChange[i]=0.0;
		}
	}
}
