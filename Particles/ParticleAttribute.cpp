/*
@file ParticleAttribute
@author John Hart, Ed Bachta, Wen Su
Particles is too big, it separates into ParticleAttribute and ParticleBehavior
*/


#include "Particles.h"
#include "ParticleStuff.h"
#include "ParticleAttribute.h"


static char *attr_xpm[] = {
	/* columns rows colors chars-per-pixel */
	"16 15 4 1",
	"  c None",
	". c #000000",
	"o c #808080",
	"X c #000080",
	/* pixels */
	"                ",
	" ...........    ",
	" .XXXXXXXXX..   ",
	" .Xoooooooo.X.  ",
	" .XXXXXXXXX.XX. ",
	" .Xoooooooo.....",
	" .XXXXXXXXXXXXX.",
	" .XoooooooooooX.",
	" .XXXXXXXXXXXXX.",
	" .XoooooooooooX.",
	" .XXXXXXXXXXXXX.",
	" .XoooooooooooX.",
	" .XXXXXXXXXXXXX.",
	" ...............",
	"                "
};

// Default contructor
ParticleAttribute::ParticleAttribute(Particles *ps, const  std::string &name)
: ParticleStuff(ps,name)
{
	
	xpm = attr_xpm;
//	perParticle = false; // sets the default value.  Child class constructors should overwrite this with the correct value;
}

// Removes current name from attributes index list
// changes the attribute's name and readds it into the list.
// The addMeTo method automatically guarantees the new name is unique.
void ParticleAttribute::setName(const std::string &new_name) {
	ParticleAttributes::iterator it = ps->attributes.find(name);
	ps->attributes.erase(it);
	name = new_name;
	addMeTo(ps);
}

void ParticleAttribute::addMeTo(Particles *new_ps)
{
	// So help me, if new_ps is null, I'm bringing the app down.
	assert(new_ps);

	std::string basename = name;
	char counter = '1';
	if (name[name.length()-3] == '(' && name[name.length()-1] == ')') {
		counter = name[name.length()-2];
		basename = name.substr(0,name.length()-3);
	}

	// if the attribute already exists and it isn't us, we need to rename ourself
	ParticleAttribute *exists = new_ps->attributes[name];
	while (exists && exists != this)
	{
		name = basename + "(" + ++counter + ")";
		exists = new_ps->attributes[name];
	}
	new_ps->attributes[name] = this;

	// Warning! Can't attach attributes here because when loading it could create things that might later collide.

	// attach all attributes needed
	//attachAttributes();

	// If the particle system is already populated, and the new attribute holds per-particle values, then we
	// need to call particleAdded the get that attribute to 'catch up' to the number of particles in the system.
	// Check to see if this new attribute holds values per-particle or just one attribute per particle system.
	// 
	// should I check if this attribute already has members?  Ifit does, how did they get there?
//	int szz = perparticle.size();
//	int szzz = new_ps->size();
 	if( perparticle.size() ) // this, I think, just checks if there are per-particle papameters in this objects ParticleStuff parent object.
	{
		// Make sure that the new attribute has allocated enough per-particle space.		
		unsigned int sz = new_ps->size();
		for(unsigned int x = 0; x < sz; x++)
			this->particleAdded();
	}

}

void ParticleAttribute::removeMe()
{
	assert(ps);
	ps->attributes.erase(ps->attributes.find(name));
}

void ParticleAttribute::attachTo(ParticleStuff *stuff)
{
//	AttachedAttribute *aa = new AttachedAttribute(stuff,this);
//	aa->shortname = aa->name = getName();
}
