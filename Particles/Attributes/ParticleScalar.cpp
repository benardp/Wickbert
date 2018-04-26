/**
 * Implementation of ParticleScalar.
 * @file ParticleScalar.cpp
 * @author Tony Kaap
 */

#include "ParticleScalar.h"

REGISTER_PARTICLESTUFF(ParticleScalar,"Attribute:ParticleScalar");

/**
 * Add particle scalar to a system of particles.
 * @param ps   The owning particle system.
 * @param name The name of this object.
 */
ParticleScalar::ParticleScalar(Particles *ps, const std::string& name, const bool& printFlag)
	: ParticleAttribute(ps, name)
{
// currently, this attribute never gets initialized to the particle system.
	// The vector t never gets allocated for the space of each particle in the system,
	// so accessing it will cause a memory error and crash wickbert.
	// At creation time, ps may be passes as NULL.  I don't know if
	// there is an `attachment` time that I can use to initialize the attribute

	if (ps) // currently, ps is NULL when this constructor is run, so there is no opportunity to initialize t
		t.resize(ps->size()); // this never gets called when the attribute is loaded after the particle system is populated

	// adding this line would mean that every particleScalar attribute will be printed in a save file
	new PSParamDoublePerParticle(this,&t,"scalarPP","per particle scalar","A per particle scalar valued attribute", printFlag);

	new PSParamInt(this,&numElements,0,"numElements","Number of elements", "Number of elements in this attribute");

	new PSParamString(this,&rbfFilename,"rbfDefaultValues.dat","rbfFileNm","rbf values filename","the name of the file that holds the values for rbf constraints - used only for RBFs");
	
	str = new std::string;
	*str = "hello";

	//new PSParRefParam(this,&par_ref,str,"<emptyTK>","PSParRefParam","Particle Reference","This is a test parameter to allow particleStuff objects to link to Particles objects, even if they are not in this object's own Particle system");

	changed = true;

	//new PSParamDoublePerParticle(this,&t,"age","particle age","Current lifetime of particle");
	//this->reset();

	//perParticle = true;
	// Although this object has been created, `ps` is always NULL when it's passed in, so t is never resized properly.
	// This flag gets set so that elements using this attribute will know that it is not currently ready for use.
	//ready = false;
}

void ParticleScalar::loadFromFile(std::string fileName)
{
 int test = 0;
}

void ParticleScalar::PSPtest()
{
	if(par_ref)
		std::string nm = par_ref->getName();

}

void ParticleScalar::reset()
{
	for(unsigned int i=0; i<t.size(); i++)
		t[i] = 0;
	changed = true;
}

void ParticleScalar::clear()
{
	t.clear();
	changed = true;
}

/**
 * Add a scalar corresponding to the new particle.
 * @param i Index of the new particle.
 */
void ParticleScalar::particleAdded()
{
	t.push_back(0.0);
	numElements = (int)t.size();
	changed = true;
}

/**
 * Callback for particle removal.
 * @param i Index of particle to be removed.
 * @see Particles::particleRemoved
 */
void ParticleScalar::particleRemoved(unsigned int i)
{
	t[i] = t.back();
	t.pop_back();
	changed = true;
}

std::vector<double> ParticleScalar::getScalar()
{
	return t;
}
double ParticleScalar::getScalar(int i)
{
//	PSPtest();
	return t[i];
}

void ParticleScalar::setScalar(const int& i, const double& value)
{
int j = t.size();
	t[i] = value;
	changed = true;
}
