/**
 * Declaration of an attribute that describes the particle connectivity on a
 * contour
 * @file Contours.h
 * @date 07/17/2006
 * @author Matei N. Stroila
 * @remarks
 */


#ifndef CONTOURS_H
#define CONTOURS_H

#include "ParticleAttribute.h"

typedef std::vector<unsigned int> Loop;
typedef struct neighbor
{
	// the first index is implicitly the index
	// the second index is stored here.
	// there are two neighbors
	unsigned int n1,n2;
} Neighbors;

class Contours : public ParticleAttribute {
	
public:
	
	MAKE_PARTICLESTUFF_NAME();
	
	Contours(Particles *ps=NULL, const std::string& name=std::string("Contours"));

	//The contours are a collection of loops. Each loop consist in the 
	//indices of the particles making the loop
	std::vector<Loop> loops;
	
	//for each particle store its neighbors
	std::vector<Neighbors> chains;
	
	//find the loop l containing the particle i
	void getLoop(const unsigned int i, Loop &l);

	//find the index j of the loop l containing the particle i
	void getLoop(const unsigned int i, unsigned int &j);

	//find the neighbors of the particle i
	void getNeighbors(const unsigned int i , Neighbors &nbs){
		nbs = chains[i];
	}
	
	void setChains(unsigned int i, const Neighbors &nbs)
	{
		chains[i]=nbs;
	}
	
	/// Callback for particle addition.
	void virtual particleAdded(); 
	
	/// Callback for particle removal.
	void virtual particleRemoved(unsigned int i);
	
	
private:
	
	Loop::iterator loop_iter;
	std::vector<Loop>::iterator loops_iter;
	
};

#endif
