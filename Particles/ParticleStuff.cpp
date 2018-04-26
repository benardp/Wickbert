

#include "Particles.h"
#include "ParticleAttribute.h"
#include "ParticleStuff.h"

const std::string ParticleStuff::getClass() { return registry_name; }


void ParticleStuff::setName(const std::string &new_name) {
	std::string old_name = name;
	name = new_name;

#if 0
	/* If this is an attribute, need to change the map */
	ParticleAttribute *me_attr = dynamic_cast<ParticleAttribute *>(this);
#if 0
	if (me_attr) {
		ParticleAttributes::iterator it;
		for (it = this->ps->attributes.begin(); it != this->ps->attributes.end(); it++) {
			if (it->second == this) {
				it->first = new_name;
			}
		}
	}
#else
	if (me_attr) {
		ParticleAttributes::iterator it = ps->attributes.find(old_name);
		ps->attributes.erase(it);
		ps->attributes[new_name] = me_attr;
	}
#endif
#endif
}

std::string ParticleStuff::getName() {return name;}

static char *stuff_xpm[] = {
	/* columns rows colors chars-per-pixel */
	"16 15 3 1",
	"  c None",
	". c #000000",
	"X c #808080",
	/* pixels */
	"                ",
	" ...........    ",
	" .XXXXXXXXX..   ",
	" .XXXXXXXXX.X.  ",
	" .XXXXXXXXX.XX. ",
	" .XXXXXXXXX.....",
	" .XXXXXXXXXXXXX.",
	" .XXXXXXXXXXXXX.",
	" .XXXXXXXXXXXXX.",
	" .XXXXXXXXXXXXX.",
	" .XXXXXXXXXXXXX.",
	" .XXXXXXXXXXXXX.",
	" .XXXXXXXXXXXXX.",
	" ...............",
	"                "
};

ParticleStuff::ParticleStuff(Particles *p, const std::string &new_name)
:ps(p), name(new_name)
{
	xpm = stuff_xpm;
}

void ParticleStuff::printContent() { std::cout << name << std::endl; }

void *
ParticleStuff::getAttributeVariable(std::string ref)
{
	unsigned int colon = ref.find(':');
	if (colon == ref.length())
		return NULL;
	std::string attr_name = ref.substr(0,colon);
	std::string param_name = ref.substr(colon+1);
	ParticleAttribute *attr = ps->getAttributeGeneric(attr_name);
	if (!attr) return NULL;
	return attr->getVariable(param_name);
}

void *
ParticleStuff::getVariable(std::string ref)
{
	for (unsigned int i = 0; i < params.size(); i++) {
		if (params[i]->shortname == ref) {
			return params[i]->ref();
		}
	}

	for (unsigned int i = 0; i < perparticle.size(); i++) {
		if (perparticle[i]->shortname == ref) {
			return perparticle[i]->ref();
		}
	}

	return NULL;
}

int ParticleStuff::qlen() {return 0;}

void ParticleStuff::getq(double *) {}

void ParticleStuff::setq(double *) {}

void ParticleStuff::qname(char **) {}

char *ParticleStuff::qtip(int i) {return "";}

char *ParticleStuff::tip() {return "";}

int ParticleStuff::qlenpp() {return 0;}

void ParticleStuff::getqpp(double *, int) {}

void ParticleStuff::setqpp(double *, int) {}

void ParticleStuff::qnamepp(char **) {}

/**
setParticleSystem should be called at constructor.
A constructor first set the default parameters.
Then call addMeTo and call attach attributes for particle behavior and shaders.
This will attach the linked attributes like implicit interrogator.
*/

void ParticleStuff::setParticleSystem(Particles *new_ps)
{
	if (!new_ps)
		return;
	ps = new_ps;
	addMeTo(ps);
}

void ParticleStuff::qname(NameVector &names)
{
	names.clear();
	
	int theQlen = qlen();
	names.reserve(theQlen);
	char** n = new char*[theQlen];
	
	qname(n);
	
	for(int i=0;i<theQlen;i++)
	{		
		names.push_back(n[i]);		
	}

	// should delete!
	delete[] n;
}

void ParticleStuff::qnamepp(NameVector &names)
{
	names.clear();
	
	int theQlen = qlenpp();
	names.reserve(theQlen);
	char** n = new char*[theQlen];
	
	qnamepp(n);
	
	for(int i=0;i<theQlen;i++)
		names.push_back(n[i]);		

	// should delete!
	delete[] n;
}

void ParticleStuff::getq(DoubleVector &q)
{	
	unsigned int theQlen = qlen();
	q.reserve(theQlen);
	if (q.size()!=theQlen)
		q.resize(theQlen);
	
	getq(&q[0]);
	
} // end getq(DoubleVector&)

void ParticleStuff::getqpp(DoubleVector &q, int i)
{	
	unsigned int theQlen = qlenpp();
	q.reserve(theQlen);
	if (q.size()!=theQlen)
		q.resize(theQlen);
	
	getqpp(&q[0],i);
	
}

void ParticleStuff::setq(DoubleVector &q)
{	
	if (q.size() != qlen())
		q.resize(qlen(), 0.0);
	
	setq(&q[0]);
	
} // end setq(DoubleVector&)

void ParticleStuff::error(std::string e)
{
	ps->error("[" + name + "] " + e);
}

void ParticleStuff::setqpp(DoubleVector &q, int i)
{	
	if (q.size() != qlenpp())
		q.resize(qlen(), 0.0);
	
	setqpp(&q[0],i);
	
}

void ParticleStuff::qshortname(NameVector &q)
{
	for(int i=0;i<qlen();++i)
	{
		std::ostringstream number;
		number << "shortname" << i;
		q.push_back(number.str());
	}
}

/** Outputs a particlestuff of the form "classname [name] [ attributes ] { <parameters> } \n".
 */
 std::ostream &operator<<(std::ostream &out, ParticleStuff *stuff) {
	std::string fullclassname = stuff->getClass();
	std::string classname = fullclassname.substr(fullclassname.find(":")+1);
	 
	out << "\t\t" << fullclassname;

	if (classname != stuff->name)
		out << " " << '"' << stuff->name << '"';

	out << " " << stuff->attachedattributes;

	out << " " << stuff->params;

	// tk[6.22.05] this outputs the non per-particle paramters of a particle system,
	// I want to implement: 
	//out << " " << stuff->perparticle;
	// Nope. We don't want to save per-particle info. If necessary, we'll record
	// it in a surface mesh and save it that way. -jch[12.27.2005]

	out << std::endl;

	return out;
}
