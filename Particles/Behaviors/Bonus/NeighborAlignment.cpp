#include "NeighborAlignment.h"


REGISTER_PARTICLESTUFF(NeighborAlignment,"Behavior:NeighborAlignment");

// Default to the X axis
NeighborAlignment::NeighborAlignment(Particles *ps)
: ParticleBehavior(ps, std::string("NeighborAlignment"))
{
	axis = 1;
	PSParamInt(this,&axis,1,"axis","alignment axis",0,2,"1 = X, 2 = Y, 3 = Z");
}

int NeighborAlignment::qlen() 
{
	return 1;	
}

void NeighborAlignment::getq(double *q)
{
	q[0] = axis;
}

void NeighborAlignment::setq(double *q)
{
	axis = q[0];
}

void NeighborAlignment::qname(char **qn)
{
	qn[0] = "Axis";
}

// Note that attachAttribute creates an attribute if it doesn't exist
void NeighborAlignment::attachAttributes()
{
	ParticleBehavior::attachAttributes();
	attachAttribute(position,std::string("ParticlePosition"));
	attachAttribute(p_orient,std::string("ParticleOrientation"));
	attachAttribute(imp, std::string("ImplicitInterrogator"));
}

// This function could, possibly, be worse.  
// But it would be hard.
//
// This is *really* just a proof of concept.
void  NeighborAlignment::applyConstraint()
{

	unsigned int n1, n2, n3, n4, n5;
	double minDistance, curDistance;

	// Iterate over each of the particles
	for (int i = 0; i < (int)ps->size(); i++) {

		minDistance=gmGOOGOL;

		// Find the closest
		for (int j = 0; j < (int)ps->size(); j++) {

			// not the same particle
			if (i!=j)
			{
				curDistance=distance(position->getPosition(j),position->getPosition(i));
				// closest
				if (curDistance<minDistance)
				{
					n1 = j;
					minDistance=curDistance;
				}
			}
		}

		
		minDistance=gmGOOGOL;

		// Find the second closest		
		for (int j = 0; j < (int)ps->size(); j++) {

			// not the same particle
			if ((i!=j) && (n1 != j))
			{
				curDistance=distance(position->getPosition(j),position->getPosition(i));
				// closest
				if (curDistance<minDistance)
				{
					n2 = j;
					minDistance=curDistance;
				}
			}
		}

		
		minDistance=gmGOOGOL;

		// Find the third closest		
		for (int j = 0; j < (int)ps->size(); j++) {

			// not the same particle
			if ((i!=j) && (n1 != j) && (n2 != j))
			{
				curDistance=distance(position->getPosition(j),position->getPosition(i));
				// closest
				if (curDistance<minDistance)
				{
					n3 = j;
					minDistance=curDistance;
				}
			}
		}

		
		minDistance=gmGOOGOL;

		// Find the fourth closest		
		for (int j = 0; j < (int)ps->size(); j++) {

			// not the same particle
			if ((i!=j) && (n1 != j) && (n2 != j) && (n3 != j))
			{
				curDistance=distance(position->getPosition(j),position->getPosition(i));
				// closest
				if (curDistance<minDistance)
				{
					n4 = j;
					minDistance=curDistance;
				}
			}
		}

		
		minDistance=gmGOOGOL;

		// Find the fith closest		
		for (int j = 0; j < (int)ps->size(); j++) {

			// not the same particle
			if ((i!=j) && (n1 != j) && (n2 != j) && (n3 != j) && (n4 != j))
			{
				curDistance=distance(position->getPosition(j),position->getPosition(i));
				// closest
				if (curDistance<minDistance)
				{
					n5 = j;
					minDistance=curDistance;
				}
			}
		}

		// Compute the average axis of the surrounding particles
		gmVector3 axisAvg = p_orient->getAxis(n1, axis) + p_orient->getAxis(n2, axis) + 
			p_orient->getAxis(n3, axis)  + p_orient->getAxis(n4, axis)  + p_orient->getAxis(n5, axis);
		axisAvg *= 1.0 / 5.0;

		axisAvg.normalize();

		// Orient the particle by setting its axis to the average of its surroundings
		p_orient->setOrientation(i, vec_to_vec_quat(p_orient->getAxis(i, axis), axisAvg) * p_orient->getOrientation(i));
	}
}


