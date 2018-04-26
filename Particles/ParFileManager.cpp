/**
* Implementation of the particle file manager.
* @file ParFileManager.cpp
* @date July 11, 2001
* @author John C. Hart
* @note Based heavily on ImpFileManager, and both should probably be a
* specialization of some FileManager master class someday.
*/

#include "ParFileManager.h"
#include "ParticleAttribute.h"
#include "ParticleBehavior.h"
#include "ParticleShader.h"
#include "ParticleOrientation.h"

#include <vector>
#include <algorithm> // for find(first,last,entry) procedure

ParFileManager::ParFileManager()
: suffix(".par"), version("v0.2"), verbose(false)
{}

void ParFileManager::setVerbose(int v)
{
	verbose = v;
}

/**
*	Read a particle system from a file.
*	The format for an object is: "stuffclassname [name] {param = value, ...} \n"
*	\param		file	File stream to read from.
*	\param		p	Pointer to new Particles object.
*	\param		name	Name of particle system
*	\returns	True if successful.
*/
bool ParFileManager::readParticles(std::ifstream &file, Particles *&p, const std::string& name)
{	
	std::string stuffclassname, qname;
	
	pfmstatus(2,"Creating Particles: " + name);

	// Create a new Particles object and set its name to name
	if (!(p = new Particles(NULL,name)))
		return pfmerror("Error: Could not create Particles: " + name);

	// group particlestuff with {
	file >> stuffclassname;
	if (stuffclassname != "{")
		return pfmerror("Read " + stuffclassname + " but expected left brace.");

	// Read particlestuff objects up to right brace
	file >> stuffclassname;
	while (!file.eof() && stuffclassname != "}")
	{
		pfmstatus(2,"New stuff: " + stuffclassname);

		// create the particle stuff object
		std::auto_ptr<ParticleStuff> auto_stuff = NEW_PARTICLESTUFF(stuffclassname);
		// avoid the auto pointer crap for now
		ParticleStuff *stuff = auto_stuff.release();

		if (!stuff)
			return pfmerror("Expected ParticleStuff (Attribute, Behavior or Shader) but could not find: " + stuffclassname);

		// init ParticleStuff
		stuff->setParticleSystem(p);

		// get its parameters
		DoubleVector q;
		NameVector qn;
		stuff->qshortname(qn);
		stuff->getq(q);

		// group with {
		file >> qname;
		if (qname != "{")
			return pfmerror("Expected left brace at beginning of parameters for " + stuffclassname + " but got " + qname);

		// read parameters until next right brace
		file >> qname;
		while(!file.eof() && qname != "}")
		{
			NameVector::iterator qni = std::find(qn.begin(),qn.end(),qname);
			int i = qni - qn.begin();
			if (i < (int) qn.size()) {
				std::string equals;
				file >> equals;
				if (equals != "=")
					return pfmerror("Expected '=' after parameter " + qname + " in " + stuffclassname);

				double value;
				file >> value;
				q[i]=value;
			} else {
				return pfmerror("Unknown parameter " + qname + " in " + stuffclassname);
			}

			// Might be a , or a }
			file >> qname;
			if (qname == "}")
				break;
			if (qname != ",")
				return pfmerror("Expected ',' between parameters in " + stuffclassname);

			file >> qname;
		}

		// get q
		stuff->setq(q);

		file >> stuffclassname;
	}

	return true;	
}

/**
  * Write a particles object to a file.
  * The format for an object is: name [params] { <child>... }
  * @param file   File stream to use for writing.
  * @param p Pointer to object.
  * @returns	  True if successful.
*/
bool ParFileManager::writeParticles(std::ofstream &file, Particles *p)
{	
	int i;
	
	// Insert object name
	file << p->name << std::endl;
	file << "{" <<std::endl;
	std::string tab(" ");

	// Write attributes first
	ParticleAttributes::iterator attr;
	for (attr = p->attributes.begin(); attr != p->attributes.end(); attr++)
	{
		ParticleAttribute *a = (*attr).second;
		file << tab <<  a->getClass() << " { ";
		DoubleVector q;
		a->getq(q);
		NameVector qn;
		a->qshortname(qn);
		// Insert object parameters
		for(i = 0; i < a->qlen(); i++)
			file << qn[i] << " = " << q[i] << ", ";
		file << "}"<< std::endl;
	}

	// Write behaviors
	ParticleBehaviors::iterator beh;
	for (beh = p->behaviors.begin(); beh != p->behaviors.end(); beh++)
	{
		file << tab << (*beh)->getClass() << " { ";
		DoubleVector q;
		(*beh)->getq(q);
		NameVector qn;
		(*beh)->qshortname(qn);
		// Insert object parameters
		for(i = 0; i < (*beh)->qlen(); i++)
			file << qn[i] << " = " << q[i] << ", ";
		file << "}"<< std::endl;
	}

	// Write Shaders
	ParticleShaders::iterator shader;
	for (shader = p->shaders.begin(); shader != p->shaders.end(); shader++)
	{
		file << tab << (*shader)->getClass() << " { ";
		DoubleVector q;
		(*shader)->getq(q);
		NameVector qn;
		(*shader)->qshortname(qn);
		// Insert object parameters
		for(i = 0; i < (*shader)->qlen(); i++)
			file << qn[i] << " = " << q[i] << " , ";
		file << "}"<< std::endl;
	}
	
	file << "}" <<std::endl;
	
	return true;
	
}


/*!
	Read an implicit surface from file.
	\param		name	File name (with or without extension)
	\param		surface	Pointer to resulting surface.
	\returns	True on success, false on failure.
*/
bool ParFileManager::readFile(const std::string& name, Particles *&p)
{	
	std::ifstream file(name.c_str());
	std::string token;
	
	if (!file)
	{
		if (verbose)
			std::cerr << "Unable to open " << name << std::endl;
		return false;
	}
	
	file >> token;
	
	// Check version
	if (!token.compare(version))
	{
		if (verbose)
			std::cout << "Reading " << version 
			<< " implicit surface description file " 
			<< name << "..." << std::endl;
	}
	else
	{
		if (verbose)
			std::cerr << name << " is not a " 
			<< version << " file!!" << std::endl;
		return false;
	}
	
	// Read base implicit name
	file >> token;

	bool readFlag=readParticles(file, p, token);
	if (readFlag)
	{
		if (verbose)
			std::cout << "done." << std::endl;
	}
	else
	{
		if (verbose)
			std::cerr << "Unable to read particle from " << name << std::endl;
		return false;
	}

	return true;
	
} // end readFile

bool ParFileManager::readFile(const std::string& name, std::vector<Particles*>& p)
{
	std::ifstream file(name.c_str());
	std::string token;
	
	if (!file)
	{
		if (verbose)
			std::cerr << "Unable to open " << name << std::endl;
		return false;
	}
	
	file >> token;
	
	// Check version
	if (!token.compare(version))
	{
		if (verbose)
			std::cout << "Reading " << version 
			<< " implicit surface description file " 
			<< name << "..." << std::endl;
	}
	else
	{
		if (verbose)
			std::cerr << name << " is not a " 
			<< version << " file!!" << std::endl;
		return false;
	}

	bool still_reading = true;
	
	while (still_reading)
	{
		// Read base implicit name
		file >> token;
		// not a new name then stop
		if (!token.compare("end"))
			break;
		if (file.eof())
			break;

		Particles *tempPtr=NULL;
		// std::cout << token << std::endl;
		still_reading = readParticles(file, tempPtr, token);
		if (still_reading)
			p.push_back(tempPtr);
	}
	if (verbose)
		std::cout << "done." << std::endl;
	return true;
}


/*!
	Write an implicit surface to file.
	\param	 name	 File name (with or without extension)
	\param	 surface Pointer to surface.
	\returns True on success, false on failure.
*/
bool ParFileManager::writeFile(const std::string& name, Particles* &p)
{	
	std::ofstream file(name.c_str());
	
	if (!file)
	{
		if (verbose)
			std::cerr << "Unable to open " << name << std::endl;
		return false;
	}
	
	file << version << std::endl;
	
	if (verbose)
		std::cout << "Writing to " << name << "... ";
	
	bool success = writeParticles(file, p);
	
	if (success)
	{
		if (verbose)
			std::cout << "done." << std::endl;
	}
	else
	{
		if (verbose)
			std::cerr << "Unable to write surface to " << name << "!" << std::endl;
	}
	
	return success;
	
} // end writeFile

/*!
	Writes a set of implicits (a scene) to file.
	\param		name		File name (with or without extension)
	\param		surfaces	Vector of pointers to surfaces.
	\returns	True on success, false on failure.
*/
bool ParFileManager::writeFile(const std::string& name, std::vector<Particles *>& p)
{
	std::ofstream file(name.c_str());
	
	if (!file)
	{
		if (verbose)
			std::cerr << "Unable to open " << name << std::endl;
		return false;
	}
	
	file << version << std::endl;
	
	if (verbose)
		std::cout << "Writing to " << name << "... ";

	bool success = true;

	for (unsigned int i = 0; i < p.size(); i++)
	{
		success = success && writeParticles(file, p[i]);
	}

	// end mark
	file << "end" << std::endl;

	if (success)
	{
		if (verbose)
			std::cout << "done." << std::endl;
	}
	else
	{
		if (verbose)
			std::cerr << "Unable to write surface to " << name << "!" << std::endl;
	}
	
	return success;
}

bool ParFileManager::writePerParticle(const std::string& name, Particles *p)
{
	std::ofstream file(name.c_str());
	
	if (!file)
	{
		if (verbose)
			std::cerr << "Unable to open " << name << std::endl;
		return false;
	}
	
	file << "# generated by particle system " << version << std::endl;
	
	if (verbose)
		std::cout << "Writing to " << name << "... ";
	
	bool success = writePerParticle(file, p);
	
	if (success)
	{
		if (verbose)
			std::cout << "done." << std::endl;
	}
	else
	{
		if (verbose)
			std::cerr << "Unable to write surface to " << name << "!" << std::endl;
	}
	
	return success;
}

bool ParFileManager::writePerParticle(std::ofstream &file, Particles *p)
{
#if 0
	file << "# particle name: " << p->name <<std::endl;
	// get normal
	ParticleOrientation *orient=p->getAttribute<ParticleOrientation>(std::string("ParticleOrientation"));
	ParticlePosition *position=p->getAttribute<ParticlePosition>(std::string("ParticlePosition"));
	if (orient)
	{
		for(unsigned int i=0;i<p->size();i++)
		{
			gmVector3 nrml = orient->getNormal(i);
			file << "x " << position->getPosition(i)[0]<< " " << position->getPosition(i)[1]<< " " << position->getPosition(i)[2]<<std::endl;
			file << "xn " << nrml[0]<< " " << nrml[1]<< " " << nrml[2]<<std::endl;
		}
	}
	else
	{
		for(unsigned int i=0;i<p->size();i++)
		{
			file << "x " << position->getPosition(i)[0]<< " " << position->getPosition(i)[1]<< " " << position->getPosition(i)[2]<< std::endl;
		}
	}
#endif
	return true;
}

