
#include "Particles.h"
#include "ParticleStuff.h"
#include "ParticleSystem.h"
#include "ParticleScalar.h"
#include "OBJPosition.h"
#include "ParticleStuffParameters.h"

#include <fstream>
#include <iostream>
#include <boost/lexical_cast.hpp>


PSParam::PSParam(ParticleStuff *parent,const std::string& sn,const std::string& n,const std::string& d = "")
{
	stuff = parent;
	shortname = sn;
	name = n;
	desc = d;
	bounded = false;
	stuff->params.push_back(this);
}

PSParamPerParticle::PSParamPerParticle(ParticleStuff *parent,const std::string& sn,const std::string& n,const std::string& d = "", const bool& pFlag)
{
	printFlag = pFlag; // Should this paramters be written out in a save file? Defaults to false
	stuff = parent;
	shortname = sn;
	name = n;
	desc = d;
	stuff->perparticle.push_back(this);
}


void
PSParamString::set(std::string s)
{
	*pv = s;

	std::string name = stuff->getName();
	std::string myParamName = this->shortname;
	
	if(myParamName == "objFileNm") // this code is a small hack to read in vertex position points from a file
	{
		OBJPosition *pPosition = stuff->ps->getAttribute<OBJPosition>(name);
		if(pPosition)
		{
			pPosition->changed = true;

//			int xVecSize = pPosition->x.size();
			pPosition->x.clear();

			char buffer[256];
			std::string path = "./../samples/"; // this will open files in the samples directory

			path.append(s);
			std::ifstream inStream (path.c_str());
			
			bool vFinished = false;
			bool vStarted = false;
			//int vertsAdded = 0;
			if(inStream.is_open())
			{
				int count = 0;
				while( !(vStarted && vFinished))
				{
					inStream.getline(buffer,100); // junk message at beginning
					
					std::string input = buffer;
					
					if(input.substr(0,2) == "v "){
						vStarted = true;
						count++;
						//if(vertsAdded < xVecSize)
						//{
							input.erase(0,2); // strip off leading characters
							double vX = 5, vY = 5, vZ = 5;
							int ws = input.find_first_of(" ",0);
							std::string xStr, yStr, zStr;
							xStr = (input.substr(0,ws));
							int ws2 = input.find_last_of(" ");
							yStr =  (input.substr(ws+1,ws2-ws-1)) ;
							vX = atof( xStr.c_str() );
							//input.erase(0,ws+1);
							vY = atof(yStr.c_str() );
							//std::string	input2 = input;
							//input2.erase(0,2);
							//input.erase(0,ws2+1);
							
							//std::string a = "0123456789";
							//a.erase(0,2);
							//a.erase(0,2);
							//int ws3 = input.find_first_of(" ",0);
							zStr = (input.substr(ws2,input.size()-ws2));
							vZ = atof( zStr.c_str() );
							gmVector3 vec(vX,vY,vZ);
							stuff->ps->addParticle();
							pPosition->setPosition(pPosition->x.size()-1,vec);
							//vertsAdded++;
						//}
					}
					else if(vStarted)
						vFinished = true;
				
				}
				//check count

				
				int h = pPosition->getNumElements();
				for(int m = 0; m < h; m++)
				{	
					gmVector3 qqq = pPosition->getPosition(m);
//					double ax = qqq[0];
//					double ay = qqq[1];
//					double az = qqq[2];
				}
				
			} //if inStream.is_open
		} // if pPosition
	} // if objFileNm

	if(myParamName == "rbfFileNm") // this code is a small hack to read in scalar values from a file
	{
		ParticleScalar * pScalar = stuff->ps->getAttribute<ParticleScalar>(name);
		if(pScalar)
		{
			pScalar->changed = true;
		
			// update the value vector with the first NUM values from the file, where NUM is the current number of particles
			//std::ifstream * inFile = new std::ifstream();
			char buffer[256];
			//std::string fileExt = ".txt";
			std::string path = "./../samples/"; // this will open files in the samples directory
			
			//add code to check if s is an output filename
			if( s.substr(0,3) != "OUT") // input file
			{
				path.append(s);
				//fileName.append(fileExt);
				std::ifstream inStream (path.c_str());
				if (inStream.is_open())
				{
					// just catch the first line
					//while(! (inStream.eof() ))
					//{
						inStream.getline(buffer,100);
						//std::string line = buffer;
						int numElementsInFile = atoi(buffer);
						int vectorSize = pScalar->getNumElements();
						int max;
						if(numElementsInFile <= vectorSize)
							max = numElementsInFile;
						else
							max = vectorSize;

						for(int i = 0; i < numElementsInFile; i++)
						{
							inStream.getline(buffer,100);
							//line = buffer;
							if(i < vectorSize)
								pScalar->setScalar(i, atof(buffer));
							else
							{
								pScalar->particleAdded();
								pScalar->setScalar(i,atof(buffer));
								vectorSize++;
							}
						}
					inStream.close();
				}
				else
					//fail - file not openeed
					;
//					int i = 0;// dummy line
					
			}else //output file - file names beginning with 'OUT' will not read in, they will write the current values of the particleScalar
				// to the file without OUT at the beginning of the name.  I.E.  s = "OUTtemp.txt" - this will output to "temp.txt"
			{
				std::string tmp;
				tmp =  s.substr(3,s.length()-3);
				path.append(tmp);
				std::ofstream outStream(path.c_str());
				if(outStream)
				{
					int numValues = pScalar->getNumElements();
					outStream << numValues << std::endl;
					for(int j = 0; j < numValues; j++)
					{
						outStream << pScalar->getScalar(j);
						if(j+1 < numValues) // prevents end of file newline
							outStream << std::endl;
					}
					outStream.close();
				}
			}
			////inFile->open(s.c_str(), ios::in);
			//if(inFile)
			//{
			//	std::string temp >> inFile;
			//	inFile->close();
			//}
			//
		} // if pScalar
	}// if PsParam name

	stuff->attachAttributes();
}

/*
PSParamComboBox::PSParamComboBox(	ParticleStuff *parent, Callback *callback, 
										const std::string& sn, 
										const std::string& n, 
										const std::string& d
									) 
									: PSParam(parent,sn,n,d)
{ 	
	//parent->params.push_back(this);
	//this line leads to bugs, because
	//here the push operation is done done by PSParam... 
	//It is NOT the same behavior as in SURFACE PARAMS!!!!!!!
	//At long term this should be unified. 
	//The solution here is in my opinion a little nicer...
	//but you have to know it! - ELMAR

	func = callback;
};
*/
PSParamComboBox::PSParamComboBox(	ParticleStuff *parent, Callback *callback, 
									const std::string& sn, 
									const std::string& n, 
									const std::string& d
								) 
								: PSParam(parent,sn,n,d)
{ 	
	func = callback;
};

PSParamComboBox::Callback::Callback(){
	//we init selected with 0 value;
	_selected=0;
};

unsigned int PSParamComboBox::Callback::nbChoices() const
{
	return _choices.size();
}

const std::string & PSParamComboBox::Callback::getChoiceString(unsigned int i) const
{
	assert(_choices.size()>0);
	if (i<_choices.size())
		return _choices[i];
	else
	{
		std::cerr<<"COMBOBOX: ERROR out of bounds for getchoicestring "<<i<<'\n';
		return _choices[0];
	}
}

const void PSParamComboBox::Callback::itemselectFromInterface(const std::string & selection)
{
	assert(_choices.size()>0);
	unsigned int i;
	for (i=0;i<_choices.size();++i)
	{
		if (_choices[i]==selection)
			break;
	}
	if (i<_choices.size())
	{
		_selected=i;
		itemselected(i);
	}
	else
	{
		std::cerr<<"COMBOBOX: InterfaceSelectedSomething not in choices! THIS SHOULD NEVER HAPPEN!"<<'\n';
		std::cerr<<"Either it was not the interface calling this funciton or it is a bug in SurfParamCombo! Selection:"<<selection<<'\n';
		std::cerr<<"choices has "<<_choices.size()<<" elements, that are:"<<'\n';
		for (unsigned int i=0;i<_choices.size();++i)
			std::cerr<<_choices[i]<<'\n';
	}
}

void PSParamComboBox::Callback::setSelectedItem(unsigned int i)
{
	if (i<_choices.size())
	{
		//accept selection and call the itemselect function
		_selected=i;
		itemselected(_selected);
	}
	else
		std::cerr<<"COMBOBOX: called item "<<i<<" which is outside choices"<<'\n';

}
unsigned int PSParamComboBox::Callback::getSelected() const
{	
	return _selected;
}

const std::string & PSParamComboBox::Callback::getSelectedString() const
{	
	return _choices[_selected];
}

std::string PSParamComboBox::get()
{

	if (func) 
	{
		//char buf[256]; sprintf(buf,"%u",func->getSelected()); return std::string(buf); 
		return std::string("")+=boost::lexical_cast<std::string>(func->getSelected());
	}
	else
		return "0";
}	

void PSParamComboBox::set(std::string itemNumber)
{
	if (!func)
		return;

	int index(atoi(itemNumber.c_str()));
	unsigned int uindex=(unsigned int)index;
	func->setSelectedItem(uindex);
}

std::ostream & PSParamComboBox::put(std::ostream &out)
{
	return out << shortname << " = " << get();
}
PSParam & PSParamComboBox::operator<<(std::string val)
{
	set(val);return (*this);
}

void PSParRefParam::set(std::string s)
{
//		std::string temp = s; // set particles to the attribute with the matching name, etc.  Run a test to see if you can pull data from the particles object
		*pv = s;
		//ImplicitInterrogator *ii = dynamic_cast<ImplicitInterrogator *>(stuff);
		//if (ii) ii->setImplicit(s);
		
		//*pv = s;
		*p = NULL;
		*p =stuff->ps->particleSystem->findParticles(s);	
		std::string nm = stuff->ps->particleSystem->name;
		if(*p) {
			std::string nmm = (*p)->name;
//			int zzz = 0; // remove this
		}
//		int removeThis = 0; //

//		ps->ParticleSystems->findParticles(s);
		//*p = FIND Particles object by name "s" 

		//p = p;
		//*p = stuff->ps->particleSystem->findParticles(*pv);
		//std::string nm;
		//if(p)
		//*p = toPSystems->findParticles(*pv);
		//stuff->attachAttributes();
}

std::string PSAttrRefParam::get()
{
	return _interfaceString;
}
void PSAttrRefParam::set(std::string s)
{
	// this test is needed!
	// The particle attributes are all initialized with a zero ps in the beginning.
	// Therefore the program would crash during this dummy step.
	if (stuff->ps==0) return;
	//init attr to null, it will be replaced
	//otherwise the user can check for success with a zero test
	*pv=0;
	
	int colon = s.rfind(':');

	// supposedly it is an attribute of the same class of particles
	//otherwise a colon would be in the name
	if (colon == -1)
	{
		*pv = stuff->ps->getAttributeGeneric(s);
	}
	else
	{
		std::string firstname = s.substr(0,colon);
		std::string secondname = s.substr(colon+1);
		colon= firstname.rfind(':');
		if (colon==-1)
		{
			//we only found one single colon, thus we are in the same ParticleSystem
			Particles * p = stuff->ps->particleSystem->findParticles(firstname);
			if (p)
			{
				*pv = p->getAttributeGeneric(secondname);
			}
		}
		else
		{
			//we have found two ':' therefore we change the particleSystem and the particles object
			//this is not very coherent, but the particle systemS access directly a particles object.
			//the syntax is the same, it separates the colon itself.
			Particles * p=stuff->ps->particleSystem->particleSystems->findParticles(firstname);
			if (p)
			{
				*pv = p->getAttributeGeneric(secondname);
			}
		}
	}
	if (*pv==0)
	{
		_interfaceString=std::string("invalid");
	}
	else
	{
		_interfaceString=s;
	}
}
PSAttrRefParam::PSAttrRefParam(	ParticleStuff *parent, ParticleAttribute **v, const std::string & init,
								const std::string& sn,const std::string& n,const std::string& d)
								: PSParamType<ParticleAttribute*>(parent,v,0,sn,n,d) 
{
	set(init); 
	defstr = _interfaceString=init;
}


//******************
//PSParRefParam::PSParRefParam(ParticleStuff *parent, Particles **par_ref, std::string s,
//						std::string sn, std::string n, std::string d)
//{
//	shortname = sn;
//	name = n;
//	desc = d;
//	p = par_ref;	// store a local reference to the surface's variable pointing to the particles object
//	*p = NULL;		// initialize the surface variable to point to a NULL particles object
//	*pv = defstr = s;
//	parent->params.push_back(this);
//}

//void PSParRefParam::attach(ParticleSystems *fromPSystems, ParticleSystems *toPSystems)
//{
//	*p = NULL;
//	p = p;
//	*p = toPSystems->findParticles(*pv);
	// [tk 5.19.05] how do I know if *p found anything or not?
	// ok, loading the pig test particles and running tests shows that attaching an rbf to the 'verts' system brings up a particle system with 6
	// particles
	//
	// that's exactly what I want, now I just have to link it properly
//}
//******************






#if 0
void
PSParamDoubleRef::set(std::string s)
{
	// check to see if constant
	char buf[256];

	if (sscanf(buf,"%g",s.c_str())) {
		*pv = buf;
		return;
	}

	// perhaps it is an attribute:variable
	std::string attr_name = s.substr(0,s.find(':'));
	std::string param_name = s.substr(s.find(':')+1);
	ParticleAttribute *attr = ps->attributes[attr_name];
	if (!attr) return;
	// found attribute now find variable
	for (int i = 0; i < attr->params.size(); i++) {
		if (attr->params[i]->shortname == param_name) {
			data = (DoubleVector *)attr->params[i]->ref();
			break;
		}
	}
	for (int i = 0; i < attr->perparticle.size(); i++) {
		if (attr->perparticle[i]->shortname == param_name) {
			data = (DoubleVector *)attr->perparticle[i]->ref();
			break;
		}
	}
}
#endif


void PSParamDoublePerParticle::set(int i,std::string s)
 {
		(*pv)[i] = atof(s.c_str()); /* Set some kind of `value changed` flag here */
		//ParticleStuff * stf = this->stuff;
		//std::string t = stf->getName();

		std::string t = this->stuff->getName();
		ParticleScalar * pScalar = this->stuff->ps->getAttribute<ParticleScalar>(t);
		if(pScalar)
			pScalar->changed = true;
	 }


void
PSParamgmVector3::set(std::string s)
{
	double x,y,z;
	sscanf(s.c_str(),"(%lf,%lf,%lf)",&x,&y,&z);
	pv->assign(x,y,z);
}

void
PSParamgmVector4::set(std::string s)
{
	double x,y,z,w;
	sscanf(s.c_str(),"(%lf,%lf,%lf,%lf)",&x,&y,&z,&w);
	pv->assign(x,y,z,w);
}

std::string
PSParamgmVector3PerParticle::get(int i)
{
	char buf[256];
	sprintf(buf,"(%g,%g,%g)",(*pv)[i][0],(*pv)[i][1],(*pv)[i][2]);
	return std::string(buf);
}

void
PSParamgmVector3PerParticle::set(int i, std::string s)
{
	double x,y,z;
	sscanf(s.c_str(),"(%lf,%lf,%lf)",&x,&y,&z);
	(*pv)[i].assign(x,y,z);
}

std::string
PSParamgmVector4PerParticle::get(int i)
{
	char buf[256];
	sprintf(buf,"(%g,%g,%g,%g)",(*pv)[i][0],(*pv)[i][1],(*pv)[i][2],(*pv)[i][3]);
	return std::string(buf);
}

void
PSParamgmVector4PerParticle::set(int i, std::string s)
{
	double x,y,z,w;
	sscanf(s.c_str(),"(%lf,%lf,%lf,%lf)",&x,&y,&z,&w);
	(*pv)[i].assign(x,y,z,w);
}

std::string
PSParamgmMatrix3PerParticle::get(int i)
{
	char buf[256];
	sprintf(buf,"((%g,%g,%g),(%g,%g,%g),(%g,%g,%g))",
		(*pv)[0][0][0],(*pv)[0][0][1],(*pv)[0][0][2],
		(*pv)[0][1][0],(*pv)[0][1][1],(*pv)[0][1][2],
		(*pv)[0][2][0],(*pv)[0][2][1],(*pv)[0][2][2]);
	return std::string(buf);
}

void
PSParamgmMatrix3PerParticle::set(int i, std::string s)
{
	double m00,m01,m02,m10,m11,m12,m20,m21,m22;
	sscanf(s.c_str(),"((%lf,%lf,%lf),(%lf,%lf,%lf),(%lf,%lf,%lf))",
		&m00,&m01,&m02,&m10,&m11,&m12,&m20,&m21,&m22);
	(*pv)[i].assign(m00,m01,m02,m10,m11,m12,m20,m21,m22);
}

/** Outputs a parameter in the form "param = value".
 * The virtual put member does the formatting, so that string values
 * are quoted.
 */
std::ostream &operator<<(std::ostream &out, PSParam *p) {
	std::string nm = p->name;
	return p->put(out);
}

/** Outputs a param list of the form "{ param = value, param = value }".
 * If the list is empty, outputs "{ }"
 */
std::ostream &operator<<(std::ostream &out, const PSParams &pvec) {
	bool output = false;

	out << "{";
//	unsigned int q = pvec.size();
	for (unsigned int i = 0; i < pvec.size(); i++) {
		std::string defstr = pvec.at(i)->defstr;
		std::string nStr = pvec.at(i)->get();
		if (pvec.at(i)->defstr != pvec.at(i)->get()) {
			if (output) out << ",";
			out << " " << pvec.at(i);
			output = true;
		}
	}
	out << " }";

	return out;
}

//*************
//************* tk 6.22.05 copying the non per-particle output operator functions to work on
// per-particle parameter classes


/** Outputs a parameter in the form "param = value".
 * The virtual put member does the formatting, so that string values
 * are quoted.
 */
std::ostream &operator<<(std::ostream &out, PSParamPerParticle *p) {
	return p->put(out);
}

/** Outputs a param list of the form "{ param = value, param = value }".
 * If the list is empty, outputs "{ }"
 */
std::ostream &operator<<(std::ostream &out, const PSParamsPerParticle &pvec) {
	bool output = false;

	int numPrintablePPParams = 0;
	for (unsigned int i = 0; i < pvec.size(); i++) {
		if( pvec.at(i)->printFlag) numPrintablePPParams++;
	}
	if(numPrintablePPParams > 0)
	{
		out << "{";
		for (unsigned int i = 0; i < pvec.size(); i++) {
			// This line currently tests each parameter in the list of parameters for this attribute to 
			// see if they are the default value - it doesn't write out default values.
			//if (pvec.at(i)->defstr != pvec.at(i)->get()) {
			//
			// For the per-particle version, this line could check the printing flag - and not print out per-particle information if the 
			// print flag is not set.  This also takes care of the problem of a paramter-wide default string check at this level.
			if (pvec.at(i)->printFlag) {
				if (output) out << ",";
				out << " " << pvec.at(i);
				output = true;
			}
		}
		out << "}";
		
	}
	return out;
}

/** Outputs a param list of the form "[ attr = name, attr = name ]".
 * If the list is empty, outputs nothing, not even empty square brackets.
 */
std::ostream &operator<<(std::ostream &out, const AttachedAttributes &aas) {
	bool output = false;

	for (unsigned int i = 0; i < aas.size(); i++) {
		if (aas.at(i)->attr_name_default != aas.at(i)->attr_name) {
			out << (output ? ", " : "[ ") << *aas.at(i);
			output = true;
		}
	}
	if (output) out << " ]";

	return out;
}

std::ostream &operator<<(std::ostream &out, const AttachedAttribute &aa) {
	out << aa.shortname << " = \"" << aa.attr_name << '"';

	return out;
}

AttachedAttribute *
AttachedAttributes::findaa(const std::string &sn)
{
	unsigned int i;
	for ( i = 0; i < size(); i++)
		if (sn == at(i)->shortname)
			break;
	return (i == size()) ? NULL : at(i);
}

std::istream &operator>>(std::istream &in, AttachedAttributes &aas) {
	std::string token;
	// assume the "[" has already been read
	while ((in >> token) && token != "]") {
		// token is the shortname of a parameter (followed by "= attributename")
		AttachedAttribute *aa = aas.findaa(token);

		if (!aa) {
			// attached attribute not found by name
			std::cerr << "No attached attribute called \"" << token << "\"." << std::endl;
			in >> token >> token;	// eat the "= value"
			continue;
		}

		// eat the "="
		if (!(in >> token) || token != "=") {
			// file ended unexpectedly when expecting = after param->shortname
			std::cerr << "Expected \"=\" after attached attribute " << aa->shortname << "\"." << std::endl;
			continue;
		}
		if (!(in >> token)) {
			// file ended unexpectedly when expecting a value for param->shortname
			std::cerr << "Unexpected eof after attached attribute " << aa->shortname << "\"." << std::endl;
			continue;
		}
		// Trim trailing comma if it exists
		if (token.substr(token.size()-1) == ",")
			token.resize(token.size()-1);
		// Trim surrounding quotes if they exist
		if (token.substr(0,1) == "\"" && token.substr(token.size()-1) == "\"")
			token = token.substr(1,token.size()-2);

		aa->attr_name = token;
	}

	return in;
}


//**************
//**************

PSParam *
PSParams::findparam(const std::string &sn)
{
	unsigned int i;
	for ( i = 0; i < size(); i++)
		if (sn == this->at(i)->shortname)
			break;
	return (i == size()) ? NULL : this->at(i);
}

PSParamPerParticle *
PSParamsPerParticle::findparam(const std::string &sn)
{
	unsigned int i;
	for ( i = 0; i < size(); i++)
		if (sn == this->at(i)->shortname)
			break;
	return (i == size()) ? NULL : this->at(i);
}


std::istream &operator>>(std::istream &in, PSParams &pvec) {
	std::string token;
	// assume the "{" has already been read
	while ((in >> token) && token != "}") {
		// token is the shortname of a parameter (followed by "= value")
		PSParam *param = pvec.findparam(token);

		if (!param) {
			// parameter not found
			in >> token >> token;	// eat the "= value"
			continue;
		}
		// eat the "="
		if (!(in >> token) || token != "=") {
			// file ended unexpectedly when expecting = after param->shortname
			continue;
		}
		if (!(in >> token)) {
			// file ended unexpectedly when expecting a value for param->shortname
			continue;
		}
		// Trim trailing comma if it exists
		if (token.substr(token.size()-1) == ",")
			token.resize(token.size()-1);

		*param << token;
	}

	return in;
}
