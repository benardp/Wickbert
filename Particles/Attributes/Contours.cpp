/**
* Implementation of an attribute that describes the particle connectivity on a
 * contour
 * @file Contours.cpp
 * @date 07/17/2006
 * @author Matei N. Stroila
 * @remarks
 */

#include "Contours.h"

REGISTER_PARTICLESTUFF(Contours,"Attribute:Contours");

Contours::Contours(Particles *ps, const std::string& name)
: ParticleAttribute(ps, name)
{
	if (ps)
		chains.resize(ps->size());
}

//find the loop l containing the particle i
void Contours::getLoop(const unsigned int i, Loop &l) 
{
	for(loops_iter = loops.begin(); loops_iter != loops.end(); loops_iter++){
		if(std::find(loops_iter->begin(), loops_iter->end(), i) != loops_iter->end()){
			l = *loops_iter;
			break;
		}
	}
}

//find the index j of the loop l containing the particle i
void Contours::getLoop(const unsigned int i, unsigned int &j) 
{
	unsigned int k = 0;	
	for(loops_iter = loops.begin(); loops_iter != loops.end(); loops_iter++){
		if(std::find(loops_iter->begin(), loops_iter->end(), i) != loops_iter->end()){
			j = k;
			break;
		}
		k++;
	}
}


/**
* Add a chain corresponding to the new particle.
 * @param i Index of the new particle.
 */
void Contours::particleAdded() 
{
	Neighbors nbs;
	nbs.n1 = nbs.n2 = ps->size()-1;
	chains.push_back(nbs);
}

/**
* Callback for particle removal.
 * @param i Index of particle to be removed.
 * @see Particles::particleRemoved
 */
void Contours::particleRemoved(unsigned int i) 
{
	chains[i] = chains.back();
	chains.pop_back();
}

