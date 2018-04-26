#include "Particles.h"
#include "ParticleSystem.h"
#include "ParticleShader.h"
#include "Surface/Surface.h"
#include <ctype.h> // for isspace()
#include "cleangl.h"



/// find a particles object with name and return NULL if not found
Particles *ParticleSystem::findParticles(std::string &name)
{
	for (unsigned int i=0;i<particles.size();i++)
		if (particles[i]->name==name)
			return particles[i];
	return NULL;
}

/** Outputs a ParticleSystem of the form "name {\n <particles>\n }\n".
 */
std::ostream &operator<<(std::ostream &out, ParticleSystem *ps) {
	out << '"' << ps->name << '"' << " {" << std::endl;
	for (int i = 0; i < (int)ps->particles.size(); i++)
		out << ps->particles[i];
	out << "}" << std::endl;

	return out;
}

/** Outputs the ParticleSystems of the form "ps1 ps2 ...".
 */
std::ostream &operator<<(std::ostream &out, const ParticleSystems &psystems) {
	for (unsigned int i = 0; i < psystems.size(); i++)
		out << psystems[i];

	return out;
}

std::istream &operator>>(std::istream &in, ParticleSystem *ps)
{
	std::string token;
	char c;

	while ((in >> c) && c != '}') {
		// read a quoted name into token
		// c should contain a doublequote
		if (c != '"') {
				ps->error(std::string("Read character '") + c +
					std::string("' when expecting a doublequoted Particles name."));
				return in;
		}

		// good, it's there, now put it back so
		// a newly created Particles can read it
		in.putback(c);

		Particles *p = new Particles(ps);
		in >> p;
	}

	return in;
}

std::istream &operator>>(std::istream &in, ParticleSystems *psystems)
{
	std::string token;
	char c;

	if ((in >> c) && c == '"') {

		// now read particlesystem name (including spaces) until next double quote
		token = "";
		in >> std::noskipws;	// read in spaces too
		while ((in >> c) && c != '"') token += c;
		in >> std::skipws;	// skip over whitespace from now on

		if (c != '"') {
			std::cerr << "Reached end of file while looking for closing doublequote in ParticleSystem name.";
			return in;
		}

		ParticleSystem *ps = new ParticleSystem();
		ps->name = token;
		psystems->push_back(ps);
		// read the particle system opening {
		if (!(in >> token) || token != "{") {
			ps->error("Opening bracket missing for ParticleSystem " + ps->name);
			return in;
		}
		in >> ps;
	} else {
		// read a quoted name into token
		// c should contain a doublequote
		std::cerr << "Read character '" << c << "' when expecting a doublequoted ParticleSystem name.\n";
	}

	return in;
}

std::string ParticleSystem::profile()
{
	std::string report = "ParticleSystem: " + name + "\n";
	for (unsigned int i = 0; i < particles.size(); i++) {
		report += particles[i]->profile() + "\n";
	}

	return report;
}

void *
ParticleSystem::getVariable(std::string ref)
{
	unsigned int colon = ref.find(':');
	if (colon == ref.length())
		return NULL;
	std::string particlesname = ref.substr(0,colon);
	ref = ref.substr(colon+1);

	Particles *p = findParticles(particlesname);
	if (!p) return NULL;

	return p->getVariable(ref);
}

bool
ParticleSystems::iserror() {
	if (errors.size()) return true;
	for (unsigned int i = 0; i < size(); i++)
		if (at(i)->errors.size()) return true;
	return false;
}

std::string
ParticleSystems::errorstring()
{
	unsigned int i,j;
	std::string s;

	for (i = 0; i < errors.size(); i++)
		s += *(errors[i]) + '\n';
	for (j = 0; j < size(); j++)
		for (i = 0; i < at(j)->errors.size(); i++)
			s += *(at(j)->errors[i]) + '\n';
	return s;
}

void
ParticleSystems::attachAttributes()
{
	for (unsigned int i = 0; i < size(); i++)
		at(i)->attachAttributes();
}

void
ParticleSystems::setSurfaces(Surfaces *s)
{
	for (unsigned int i = 0; i < size(); i++)
		at(i)->surfaces = s;
}

void
ParticleSystems::surfacesAttach()
{
	for (unsigned int i = 0; i < size(); i++)
		at(i)->surfaces->attach();
}

void
ParticleSystems::fullUpdate()
{
	for (unsigned int i = 0; i < size(); i++)
		at(i)->fullUpdate();
}

std::string ParticleSystems::profile()
{
	std::string report;

	for (unsigned int i = 0; i < size(); i++)
		report += at(i)->profile() + "\n";

	return report;
}

void ParticleSystems::animation()
{

	if(!animate) return;
	
	for (unsigned int j = 0; j < this->size(); j++){
		//pass the Euler angles and the zoom factor
		this->at(j)->setEulerAnglesAndZoom(xrot,yrot,zoom);
		for(unsigned int i = 0; i < this->at(j)->particles.size(); i++)
			this->at(j)->particles[i]->event(WB_CLICK);
	}
	
	if(rotateX) xrot += speed;
	if(rotateY) yrot += speed;
	if(rotateZ) zrot += speed;
	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef( 0.0f, 0.0f, (GLfloat)zoom);
	glRotatef((GLfloat)yrot,0.0f,1.0f,0.0f);
	glRotatef((GLfloat)xrot,1.0f,0.0f,0.0f);
	glRotatef((GLfloat)zrot,0.0f,0.0f,1.0f);
	
	for (unsigned int j = 0; j < this->size(); j++){
		//pass the Euler angles and the zoom factor
		this->at(j)->setEulerAnglesAndZoom(xrot,yrot,zoom);
		for(unsigned int i = 0; i < this->at(j)->particles.size(); i++)
			this->at(j)->particles[i]->event(WB_DRAG);
	}
	
	
}
