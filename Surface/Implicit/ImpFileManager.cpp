/**
* Implementation of the implicit file manager.
* @file ImpFileManager.cpp
* @date July 11, 2001
* @author Ed Bachta
*/

#include "ImpFileManager.h"

#include <set>
#include <vector>
#include "Surface/Implicit/Variational/RBF.h"


ImpFileManager::ImpFileManager() : suffix(".imp"), version("v0.3"), verbose(false)
{ }

void ImpFileManager::setVerbose(int v)
{ verbose = v; }

/**
*	Read an implicit object from file.
*	The format for an object is: name [params] { <child>... }
*	\param		file	File stream to read from.
*	\param		object	Pointer to new object.
* \param    name The name of the object being created.
*	\returns	True if successful.
*/
//bool ImpFileManager::readImplicit(std::ifstream &file, std::auto_ptr<Implicit> &object, std::string name)
bool ImpFileManager::readImplicit(std::ifstream &file, Implicit * &object, std::string &name)
{	
	std::string token;
	
	if (verbose > 1)
	{
		std::cerr << "creating " << name << std::endl;
	}

	std::auto_ptr<Implicit> temp=NEW_IMPLICIT(name);

	object = temp.release();
	
	if (object==NULL)
	{
		// Invalid object name
		if (verbose) 
			std::cerr << "Object \"" << name << "\" unknown!" << std::endl;
		return false;
	}

	//call the specialized reading function for the created implicit
	
	if(!object->readImplicit(file,verbose)) 
		return false;

/*All of this "specialized" loading activity has been moved to the individual Implicit::readImplicit virtual members -SIK
	int num_params = 0;
	file >> num_params;
	std::vector<double> params = std::vector<double>(num_params, 0.0);

	for (int i = 0; i < num_params; i++)
	{
		file >> params[i];
	}

	if (params.size() != object->qlen())
	{	// for some reason the object has more parameters
		// listed than it can accept

		if (dynamic_cast<Algebraic*>(object))
		{	// if this is an algebraic, we should increase its degree
			delete object;

			int degree = 0;
			while (Algebraic::coefficients(degree) < (int)params.size())
				degree++;

			object = new Algebraic(degree);
		}

		else if(RBF* rbf_obj = dynamic_cast<RBF*>(object))
		{
			int num_values, starting_point;
			// If this is an RBF, the rest of the values represent control point locations, normals, and
			//  function values (0 = surface, -1 = interior, +1 = exterior).

			// If it is a compact RBF, we have 2 real parameters to worry about.
			if(dynamic_cast<CompactRBF*>(object))
			{
				// Put in the real parameters.
				std::vector<double> real_params = std::vector<double>(2, 0.0);
				real_params[0] = params[0];
				real_params[1] = params[1];
				object->setq(real_params);
				
				num_values = params.size() - 2;
				starting_point = 2;
			}

			else
			{
				num_values = params.size();
				starting_point = 0;
			}

			// Pull out the rest of the parameters 7 at a time and create constraints out of them.
			//  If the number of values is not a multiple of 7, throw away whatever is left.
			if(num_values % 7)
				num_values -= (num_values % 7);

			int j = starting_point;
			rbf_obj->constraints.clear();

			while(j < params.size())
			{
				gmVector3 position, normal;
				double fun_val;
				bool normal_constraint;

				position = gmVector3(params[j], params[j + 1], params[j + 2]);
				normal = gmVector3(params[j + 3], params[j + 4], params[j + 5]);
				fun_val = params[j + 6];
				
				if(normal.lengthSquared() == 0.0)
					normal_constraint = false;
				else
					normal_constraint = true;

				rbf_obj->constraints.push_back(RBFModelerConstraint(position, normal, fun_val, normal_constraint));

				j += 7;
			}
			rbf_obj->update();

			return true;
		}

		else
		{	// dont' know what to do
			params.resize(object->qlen());
		}
	}
	
	object->setq(params);
*/

	std::getline(file, token);

	// Children included in brackets
	if (token.find('{') != std::string::npos)
	{
		std::vector<Implicit*> children;		
		file >> token; // read first child's name
		
		while(token.compare("}"))
		{			
			Implicit* child=NULL; 
			
			if (!readImplicit(file, child, token))
			{
				if (verbose)
					std::cerr << "Problem reading child of " << name << std::endl;
				return false;
			}
			
			children.push_back(child);			
			// child.release();			
			file >> token; // read next child's name			
		}
		
		// Add children to object
		if(!object->setChildren(children))
		{
			if (verbose)
				std::cerr << "Problem assigning children of " << name << std::endl;
			return false;
		}		
	}
	
	return true;
		
} // end readImplicit;


/**
 * Write an implicit object to file.
 * The format for an object is: name [params] { <child>... }
 * @param file   File stream to use for writing.
 * @param object Pointer to object.
 * @param tabs A string for the tab.
 * @returns	  True if successful.
 */
bool ImpFileManager::writeImplicit(std::ofstream &file, Implicit *&object, std::string tabs)
{	
	unsigned i=0;
	
	RBF* rbf_obj = dynamic_cast<RBF*>(object);

	// Insert object name
	file << tabs << object->name() << " ";
	std::cout << object->name() << std::endl;
	
	// Retrieve object parameters
	std::valarray<double> params(object->qlen());
	object->getq(params);

	if(rbf_obj)
	{
		int length = rbf_obj->plen() + (7 * rbf_obj->constraints.size());
		file << length << " ";
	}
	
	else
		file << object->plen() << " ";
	
	// Insert object parameters
	
	for(i=0;i<object->plen();i++)
		file << params[i] << " ";
	
	// If we have an RBF, we need to add the control point locations into the file.
	if(rbf_obj)
	{
		for(unsigned int j = 0; j < rbf_obj->constraints.size(); j++)
		{
			gmVector3 position = rbf_obj->constraints[j].c;
			gmVector3 normal = rbf_obj->constraints[j].n;
			double fun_val = rbf_obj->constraints[j].h;

			file << position[0] << " " << position[1] << " " << position[2] << " ";
			file << normal[0] << " " << normal[1] << " " << normal[2] << " ";
			file << fun_val << " ";
		}
	}

	// Retrieve children
	std::vector<Implicit*> children; 
	object->getChildren(children);
	
	if (children.size()>0)
	{		
		// Insert children
		file << '{' << std::endl;
		
		for(i=0;i<children.size();i++)
		{
			if(!writeImplicit(file, children[i], tabs+'\t'))
				return false;
		}
		
		file << tabs << '}' << std::endl;
		
	} // end child check
	
	file << std::endl;
	
	return true;
	
} // end writeImplicit;


/*!
	Read an implicit surface from file.
	\param		name	File name (with or without extension)
	\param		surface	Pointer to resulting surface.
	\returns	True on success, false on failure.
*/
bool ImpFileManager::readFile(const std::string& name, Implicit *&surface)
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

	bool readFlag=readImplicit(file, surface, token);
	if (readFlag)
	{
		if (verbose)
			std::cout << "done." << std::endl;
	}
	else
	{
		if (verbose)
			std::cerr << "Unable to read surface from " << name << std::endl;
		return false;
	}

	return true;
	
} // end readFile

bool ImpFileManager::readFile(const std::string& name, std::vector<Surface*>& surfaces)
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

		// std::auto_ptr<Implicit> tempPtr;
		Implicit *tempPtr=NULL;

		still_reading = readImplicit(file, tempPtr, token);
		if (still_reading)
			surfaces.push_back(tempPtr);
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
bool ImpFileManager::writeFile(const std::string& name, Implicit* &surface)
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
	
	bool success = writeImplicit(file, surface, std::string(""));
	
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
bool ImpFileManager::writeFile(const std::string& name, std::vector<Surface *>& surfaces)
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

	for (unsigned int i = 0; i < surfaces.size(); i++)
	{
		Implicit *imp= dynamic_cast<Implicit *> (surfaces[i]);
		if (imp)
			success = success && writeImplicit(file, imp, std::string(""));
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

