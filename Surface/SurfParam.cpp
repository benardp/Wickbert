/** @file SurfParam.cpp
 *  @author John C. Hart
 *  @date 29 Dec. 2004
 */

#include "Surface/SurfParam.h"
#include "Particles/ParticleSystem.h"
#include "Surface/Surface.h"
#include "Surface/Implicit/Implicit.h"


SurfParamBool::SurfParamBool(Surface *parent, bool *v, bool init,const std::string& sn,const std::string& n,const std::string& d)
: SurfParamType<bool>(v,init,sn,n,d) 
{ 
	parent->params.push_back(this); 
}

SurfParamButton::SurfParamButton(Surface *parent, Callback *callback, const std::string& sn, const std::string& n, const std::string& d) : SurfParam(sn,n,d) 
{ 
	parent->params.push_back(this); 
	func = callback; 
};

SurfParamComboBox::SurfParamComboBox(	Surface *parent, Callback *callback, 
										const std::string& sn, 
										const std::string& n, 
										const std::string& d
									) 
									: SurfParam(sn,n,d)
{ 	
	parent->params.push_back(this); 
	func = callback;
};


SurfParamComboBox::Callback::Callback(){
	//we init selected with 0 value;
	_selected=0;
};

unsigned int SurfParamComboBox::Callback::nbChoices() const
{
	return _choices.size();
}

const std::string & SurfParamComboBox::Callback::getChoiceString(unsigned int i) const
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

const void SurfParamComboBox::Callback::itemselectFromInterface(const std::string & selection)
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

void SurfParamComboBox::Callback::setSelectedItem(unsigned int i)
{
	if (i<_choices.size())
	{
		//accept selection and call the itemselected function
		_selected=i;
		itemselected(_selected);
	}
	else
		std::cerr<<"COMBOBOX: called item "<<i<<" which is outside choices"<<'\n';

}
unsigned int SurfParamComboBox::Callback::getSelected() const
{	
	return _selected;
}

const std::string & SurfParamComboBox::Callback::getSelectedString() const
{	
	return _choices[_selected];
}

std::string SurfParamComboBox::get()
{
	if (func) 
		return std::string("")+=(func->getSelected());
	else
		return "0";
}	

void SurfParamComboBox::set(std::string itemNumber)
{
	if (!func)
		return;

	int index(atoi(itemNumber.c_str()));
	unsigned int uindex=(unsigned int)index;
	func->setSelectedItem(uindex);
}




SurfParamInt::SurfParamInt(Surface *parent, int *v, int init,const std::string& sn,const std::string& n,const std::string& d)
: SurfParamType<int>(v,init,sn,n,d) 
{ 
	parent->params.push_back(this); 
}

SurfParamDouble::SurfParamDouble(Surface *parent, double *v, double init,const std::string& sn,const std::string& n,const std::string& d)
: SurfParamType<double>(v,init,sn,n,d) 
{ 
	parent->params.push_back(this); 
}


SurfParamString::SurfParamString(Surface *parent, std::string *v, std::string init,const std::string& sn,const std::string& n,const std::string& d)
: SurfParamType<std::string>(v,init,sn,n,d) 
{	
	parent->params.push_back(this); 
}

SurfParamgmVector3::SurfParamgmVector3(Surface *parent, gmVector3 *v, gmVector3 init, const std::string& sn, const std::string& n, const std::string& d)
: SurfParamType<gmVector3>(v,init,sn,n,d) 
{
	parent->params.push_back(this);
}


SurfStringParam::SurfStringParam(Surface *parent, std::string *s, std::string sn, std::string n, std::string d)
{
	shortname = sn;
	name = n;
	desc = d;
	v = def = *s;
	parent->params.push_back(this);
}

SurfImpRefParam::SurfImpRefParam(Surface *parent, Implicit **imp_ref)
{
	char index = 'A' + parent->params.size();
	shortname = index;
	name = std::string("Component ") + index;
	desc = "A reference by name of an Implicit.";
	imp = imp_ref;
	*imp = NULL;
	v = def = "<empty>";
	parent->params.push_back(this);
}


SurfImpRefParam::SurfImpRefParam(Surface *parent, Implicit **imp_ref, std::string s,
									std::string sn, std::string n, std::string d = "")
{
	shortname = sn;
	name = n;
	desc = d;
	imp = imp_ref;
	*imp = NULL;
	v = def = s;
	parent->params.push_back(this);
}


void
SurfImpRefParam::attach(Surfaces *surfs, ParticleSystems *psystems)
{
	*imp = dynamic_cast<Implicit *>(surfs->withName(v));
	if (*imp==0)
		v="<invalid>";
}

SurfSurfRefParam::SurfSurfRefParam(Surface *parent, Surface **surf_ref)
{
	char index = 'A' + parent->params.size();
	shortname = index;
	name = std::string("Component ") + index;
	desc = "A reference by name of an Implicit.";
	surf = surf_ref;
	*surf = NULL;
	v = def = "<empty>";
	parent->params.push_back(this);
}


SurfSurfRefParam::SurfSurfRefParam(Surface *parent, Surface **surf_ref, std::string s,
									std::string sn, std::string n, std::string d = "")
{
	shortname = sn;
	name = n;
	desc = d;
	surf = surf_ref;
	*surf = NULL;
	v = def = s;
	parent->params.push_back(this);
}


void
SurfSurfRefParam::attach(Surfaces *surfs, ParticleSystems *psystems)
{
	*surf = surfs->withName(v);
	if (*surf==0)
		v="<invalid>";
}



SurfParRefParam::SurfParRefParam(Surface *parent, Particles **par_ref, std::string s,
									std::string sn, std::string n, std::string d = "")
{
	shortname = sn;
	name = n;
	desc = d;
	p = par_ref;	// store a local reference to the surface's variable pointing to the particles object
	*p = NULL;		// initialize the surface variable to point to a NULL particles object
	v = def = s;
	parent->params.push_back(this);
}

void
SurfParRefParam::attach(Surfaces *surfs, ParticleSystems *psystems)
{
	*p = NULL;
	p = p;
	*p = psystems->findParticles(v);
	// [tk 5.19.05] how do I know if *p found anything or not?
	// ok, loading the pig test particles and running tests shows that attaching an rbf to the 'verts' system brings up a particle system with 6
	// particles
	//
	// that's exactly what I want, now I just have to link it properly
}

SurfAttrRefParam::SurfAttrRefParam(Surface *parent, ParticleAttribute **attr_ref, std::string s,
									std::string sn, std::string n, std::string d = "")
{
	shortname = sn;
	name = n;
	desc = d;
	attr = attr_ref;	// store a local reference to the surface's variable pointing to the particles attribute
	*attr = NULL;	// initialize the surface variable to point to a NULL particles attribute
	v = def = s;
	parent->params.push_back(this);
}

void SurfAttrRefParam::attach(Surfaces *surfs, ParticleSystems *psystems)
{
	int colon = v.rfind(':');

	// need at least one colon to indicate which particles holds the attribute
	if (colon == -1)
		return;

	std::string particlesname = v.substr(0,colon);
	std::string attrname = v.substr(colon+1);

	Particles *p = psystems->findParticles(particlesname);
	if (!p){
		v="<invalid>";
		return;
	}

	/* attr is a pointer to the location of a Surface member
	 * variable holding the pointer to the attribute.
	 */
	*attr = p->getAttributeGeneric(attrname);
	if (*attr==0)
		v="<invalid>";
}

/** Read in a stream of parameter settings of the form name=value.
 */
std::istream &operator>>(std::istream &in, SurfParams *s)
{
	std::string token,name,value;
	char c;

	/* so long as the next entry isn't a number (which would signify the start
	 * of a q-parameter value list) read a sequence of name=value entries, separated
	 * by commas.
	 */

	/* eat up any leading white space */
	in >> std::noskipws;
	while (isspace(in.peek()))
		in >> c;

    while (isalpha(in.peek())) {

		/* read name */
		name.clear();
		while ((in >> c) && c != '=' && !isspace(c))
			name += c;

		/* read equals sign */
		while (c != '=')
			in >> c;
		while (isspace(in.peek()))
			in >> c;

		/* read value */
		if (in.peek() == '"') {
			in >> c;
			value.clear();
			while ((in >> c) && c != '"') {
				value += c;
			}
		} else {
			in >> value;
		}

		unsigned int i;
		for (i = 0; i < s->size(); i++) {
			if (s->at(i)->shortname == name) {
				s->at(i)->set(value);
				break;
			}
		}
		if (i == s->size()) {
			std::cerr << "Parameter " << name << " not found.";
			return in;
		}

		/* eat up white space */
		while (isspace(in.peek()))
			in >> c;
	}
	
	in >> std::skipws;
	return in;
}

std::ostream &operator<<(std::ostream &out, SurfParams *s)
{
	for (unsigned int i = 0; i < s->size(); i++) {
		out << s->at(i) << "\n";
	}
	return out;
}

std::ostream &operator<<(std::ostream &out, SurfParam *param)
{
	out << "\t" << param->shortname << " = " << param->get();

	return out;
}
