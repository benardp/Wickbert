 /** @file ParticleStuffParameters.h
 *  @author John C. Hart
 */

#ifndef PSPARAM_H
#define PSPARAM_H

//ignore deprecated warning
#pragma warning(disable : 4996)

#include <string>
#include <vector>
#include <iostream> 
#include "libgm/gmVec3.h"
#include "libgm/gmVec4.h"
#include "libgm/gmMat3.h"

class ParticleStuff;
class ParticleAttribute;
class Particles;
/** Need to set up so that we can enter something like
 * AdaptiveRepulsionData:radius for a value and have it automatically
 * find the attribute and variable and load the values.
 * Will need some kind of setup() function.
 *
 */
class NamedParam
{
public:
	ParticleStuff *stuff;

	/** A one-word name for the parameter */
	std::string shortname;

	/** A one-line name for the parameter */
	std::string name;

	/** A full description of the parameter */
	std::string desc;

	virtual void *ref() { return 0; }
	
	virtual ~NamedParam(){};
};

class PSParam : public NamedParam
{
public:
	/** bounded? */
	bool bounded;

	/** default value */
	std::string defstr;

#if 0
	/** choices is a list of names and corresponding values
	 * for the user interface to provide in a dropdown list.
	 * If empty, then UI just provides a text entry.
	 */
	std::vector<std::pair<std::string,std::string>> choices;
#endif

	PSParam(ParticleStuff *parent, const std::string& sn, const std::string& n, const std::string& d);
//	{ stuff = parent; name = n; desc = d; bounded = false;
//	  stuff->param.push_back(this); }

	//output A SINGLE token! several are NOT allowed!
	virtual std::string get() = 0;
	virtual void set(std::string) = 0;

	//should it not be out<<get()???? In this way by default the value stored is in the get() function.
	//this woule make more sense...
	//but I am afraid this screws up the whole loading routine... because everything is crossing and jumping...
	// I think a work around would be to define get() as return "empty" and set as return. Then we could call it like this.
	//PUT should NOT be virtual because the loader expects a SINGLE token... this is not obvious!
	virtual std::ostream &put(std::ostream &out) { return out; }

	//NO, again by default it should be set(val); return *this;
	//otherwise we have to duplicate code everywhere!
	virtual PSParam &operator<<(std::string val) { return *this; }
	
	virtual ~PSParam(){};
};

std::ostream &operator<<(std::ostream &out, PSParam *p);

class PSParamPerParticle : public NamedParam
{
public:
	bool printFlag; // Should this paramters be written out in a save file? Defaults to false
	PSParamPerParticle(ParticleStuff *parent, const std::string& sn, const std::string& n, const std::string& d, const bool& pFlag = false);
	virtual std::ostream &put(std::ostream &out) { return out; }
	virtual std::string get(int) = 0;
	virtual void set(int, std::string) = 0;
	virtual ~PSParamPerParticle(){};
};

std::ostream &operator<<(std::ostream &out, PSParamPerParticle *p);

template<class ValType>
class PSParamType : public PSParam
{
public:
	ValType *pv;

	/** Low value of the bound */
	ValType low;

	/** High value of the bound */
	ValType high;

	PSParamType(ParticleStuff *parent, ValType *v, ValType init,const std::string& sn,const std::string& n,const std::string& d = "")
	: PSParam(parent,sn,n,d) { pv = v; *pv = init; }

	PSParamType(ParticleStuff *parent, ValType *v, ValType init, ValType l, ValType h,const std::string& sn,const std::string& n,const std::string& d = "")
	: PSParam(parent,sn,n,d) { pv = v; *pv = init; bounded = true; low = l; high = h; }

	virtual void *ref() { return pv; }

	//what about the low and high bounds ??? - Elmar
	virtual std::ostream &put(std::ostream &out) { return out << shortname << " = " << get(); }

	virtual PSParam &operator<<(std::string val) { set(val); return *this; }
	
	virtual ~PSParamType(){};
};

template<class ValType>
class PSParamTypePerParticle : public PSParamPerParticle
{
public:
	std::vector<ValType> *pv;

	PSParamTypePerParticle(ParticleStuff *parent, std::vector<ValType> *v,const std::string& sn,const std::string& n,const std::string& d = "", const bool& pFlag = false)
	: PSParamPerParticle(parent,sn,n,d,pFlag) { pv = v; }

	// I think this has to be defined here due to the template
	virtual std::ostream &put(std::ostream &out)
	{
		//bool output = false;
		int sz = (int)pv->size();
		if (sz > 0)
		{
			out << shortname << " (" << sz << ") =" << std::endl;	
			for(int i = 0; i < sz; i++)
			{
				//if (output) out << ",";
				out << "\t\t\t " << get(i) << "," << std::endl;
				//output = true;
			}
		}
		return out;
	}

	virtual void *ref() { return pv; }
	
	virtual ~PSParamTypePerParticle(){};
};

class PSParamBool : public PSParamType<bool>
{
public:
	virtual std::string get() { if (*pv) return std::string("true"); else return std::string("false"); }
	virtual void set(std::string s) { *pv = (s[0] == 'T' || s[0] == 't'); }

	PSParamBool(ParticleStuff *parent, bool *v, bool init,const std::string& sn,const std::string& n,const std::string& d = "")
	: PSParamType<bool>(parent,v,init,sn,n,d) { defstr = get(); }
	
	virtual ~PSParamBool(){};
};
/*
class PSParamBoolPerParticle : public PSParamTypePerParticle<bool>
{
public:
	virtual std::string get(int i) { char buf[32]; sprintf(buf,"%g",(*pv)[i]); return std::string(buf); }
	virtual void set(int i,std::string s) { (*pv)[i] = atof(s.c_str()); }
	
	PSParamBoolPerParticle(ParticleStuff *parent, std::vector<bool> *v,const std::string& sn,const std::string& n,const std::string& d = "")
	: PSParamTypePerParticle<bool>(parent,v,sn,n,d) {}
	
	virtual ~PSParamBoolPerParticle(){};
	
};
*/
class PSParamInt : public PSParamType<int>
{
public:
	virtual std::string get() { char buf[32]; sprintf(buf,"%d",*pv); return std::string(buf); }
	virtual void set(std::string s) { *pv = atoi(s.c_str()); }

	PSParamInt(ParticleStuff *parent, int *v, int init,const std::string& sn,const std::string& n,const std::string& d = "")
	: PSParamType<int>(parent,v,init,sn,n,d) { defstr = get(); }

	PSParamInt(ParticleStuff *parent, int *v, int init,const std::string& sn,const std::string& n, int l, int h,const std::string& d = "")
	: PSParamType<int>(parent,v,init,l,h,sn,n,d) { defstr = get(); }
	virtual ~PSParamInt(){};
};

class PSParamIntPerParticle : public PSParamTypePerParticle<int>
{
public:
	virtual std::string get(int i) { char buf[32]; sprintf(buf,"%d",(*pv)[i]); return std::string(buf); }
	virtual void set(int i,std::string s) { (*pv)[i] = atoi(s.c_str()); }
	PSParamIntPerParticle(ParticleStuff *parent, std::vector<int> *v,const std::string& sn,const std::string& n,const std::string& d = "")
	: PSParamTypePerParticle<int>(parent,v,sn,n,d) {}
	virtual ~PSParamIntPerParticle(){};
};

class PSParamDouble : public PSParamType<double>
{
public:
	virtual std::string get() { char buf[32]; sprintf(buf,"%g",*pv); return std::string(buf); }
	virtual void set(std::string s) { *pv = atof(s.c_str()); }

	PSParamDouble(ParticleStuff *parent, double *v, double init,const std::string& sn,const std::string& n,const std::string& d = "")
	: PSParamType<double>(parent,v,init,sn,n,d) { defstr = get(); }

	PSParamDouble(ParticleStuff *parent, double *v, double init, double l, double h,const std::string& sn,const std::string& n,const std::string& d = "")
	: PSParamType<double>(parent,v,init,l,h,sn,n,d) { defstr = get(); }
	
	virtual ~PSParamDouble(){};
};

#if 0
class PSParamDoubleRef : public PSParamDouble
{
public:
	std::string refstr;

	virtual std::string get() { return refstr; }
	virtual void set(std::string s);
	virtual ~PSParamDoubleRef(){};

}
#endif

class PSParamDoublePerParticle : public PSParamTypePerParticle<double>
{
public:
	virtual std::string get(int i) {char buf[32]; sprintf(buf,"%g",(*pv)[i]); return std::string(buf);}
	virtual void set(int i,std::string s);

	PSParamDoublePerParticle(ParticleStuff *parent, std::vector<double> *v,const std::string& sn,const std::string& n,const std::string& d = "", const bool& pF = false)
	: PSParamTypePerParticle<double>(parent,v,sn,n,d,pF) {}
	
	virtual ~PSParamDoublePerParticle(){};

};

//template class PSParam<std::string>;
class PSParamString : public PSParamType<std::string>
{
public:
	virtual std::string get() { return *pv; }
	virtual void set(std::string s);

	PSParamString(ParticleStuff *parent, std::string *v, std::string init,const std::string& sn,const std::string& n,const std::string& d = "")
	: PSParamType<std::string>(parent,v,init,sn,n,d) { defstr = get(); }

	PSParamString(ParticleStuff *parent, std::string *v, std::string init, std::string l, std::string h,const std::string& sn,const std::string& n,const std::string& d = "")
	: PSParamType<std::string>(parent,v,init,l,h,sn,n,d) { defstr = get(); }

	virtual std::ostream &put(std::ostream &out) { return out << shortname << " = " << '"' << get() << '"'; }

	virtual PSParam &operator<<(std::string val) { set((val[0] == '"') ? val.substr(1,val.size()-2) : val); return *this; }
	
	virtual ~PSParamString(){};

};

// Used to reference a Particles object
class PSParRefParam : public PSParamString
{
public:
	Particles **p;

	PSParRefParam(ParticleStuff *parent, Particles** par_ref, std::string *v, std::string init,const std::string& sn,const std::string& n,const std::string& d = "")
	: PSParamString(parent,v,init,sn,n,d) {p = par_ref; *p = NULL; }
//ImplicitInterrogator is a good model for this constructor

	//PSParRefParam(ParticleStuff *parent, Particles **par_ref, std::string s, std::string sn, std::string n,const std::string d= "") : PSParamString(parent, v,init,sn,n,d);

//	virtual void attach(Surfaces *surfs, ParticleSystems *psystems);
	//Surfaces is a vector container of Surface objects
	//What is a analogue for Particles? ParticleSystems? - ParticleSystems::findParticles returns Particles* - looks good
	//virtual void attach(ParticleSystems *fromPSystems, ParticleSystems *toPSystems);

	//virtual void attachAtributes(
	
	// I need to override the set() function, which gets called in Particles instead of attach();
	virtual void set(std::string s);
	
	
	virtual ~PSParRefParam() {}
};

//class creating a reference to attributes of a particle system.
//It is NOT creating these values!
//Simple example: Imagine you have a particle class, whose repulsion is not only based on its own members, but also on other particles
//Thus you want to create a reference to this other particles object -- EE
//If there is a better way, I will change it, but so far I think it is the only possibility...
class PSAttrRefParam: public PSParamType<ParticleAttribute*>
{
	std::string  _interfaceString;
public:
	//this function is overridden
	//it checks the new string and tests for attributes
	virtual void set(std::string s);
	virtual std::string get();
	
	PSAttrRefParam(ParticleStuff *parent, ParticleAttribute **v, const std::string & init,const std::string& sn,const std::string& n,const std::string& d = "");
	virtual ~PSAttrRefParam(){};
};

class PSParamgmVector3 : public PSParamType<gmVector3>
{
public:
	virtual std::string get() { char buf[1024]; sprintf(buf,"(%g,%g,%g)",(*pv)[0],(*pv)[1],(*pv)[2]); return std::string(buf); }
	virtual void set(std::string s);

	PSParamgmVector3(ParticleStuff *parent, gmVector3 *v, gmVector3 init,const std::string& sn,const std::string& n,const std::string& d = "")
	: PSParamType<gmVector3>(parent,v,init,sn,n,d) { defstr = get(); }

	PSParamgmVector3(ParticleStuff *parent, gmVector3 *v, gmVector3 init, gmVector3 l, gmVector3 h,const std::string& sn,const std::string& n,const std::string& d = "")
	: PSParamType<gmVector3>(parent,v,init,l,h,sn,n,d) { defstr = get(); }
	
	virtual ~PSParamgmVector3(){};
};

class PSParamgmVector4 : public PSParamType<gmVector4>
{
public:
	virtual std::string get() { char buf[1024]; sprintf(buf,"(%g,%g,%g,%g)",(*pv)[0],(*pv)[1],(*pv)[2],(*pv)[3]); return std::string(buf); }
	virtual void set(std::string s);

	PSParamgmVector4(ParticleStuff *parent, gmVector4 *v, gmVector4 init,const std::string& sn,const std::string& n,const std::string& d = "")
	: PSParamType<gmVector4>(parent,v,init,sn,n,d) { defstr = get(); }

	PSParamgmVector4(ParticleStuff *parent, gmVector4 *v, gmVector4 init, gmVector4 l, gmVector4 h,const std::string& sn,const std::string& n,const std::string& d = "")
	: PSParamType<gmVector4>(parent,v,init,l,h,sn,n,d) { defstr = get(); }
	
	virtual ~PSParamgmVector4(){};

};

class PSParamgmVector3PerParticle : public PSParamTypePerParticle<gmVector3>
{
public:
	virtual std::string get(int i);
	virtual void set(int i, std::string s);

	PSParamgmVector3PerParticle(ParticleStuff *parent, std::vector<gmVector3> *v,std::string sn,const std::string& n,const std::string& d = "")
	: PSParamTypePerParticle<gmVector3>(parent,v,sn,n,d) {}
	
	virtual ~PSParamgmVector3PerParticle(){};
};

class PSParamgmVector4PerParticle : public PSParamTypePerParticle<gmVector4>
{
public:
	virtual std::string get(int i);
	virtual void set(int i, std::string s);

	PSParamgmVector4PerParticle(ParticleStuff *parent, std::vector<gmVector4> *v,std::string sn,const std::string& n,const std::string& d = "")
	: PSParamTypePerParticle<gmVector4>(parent,v,sn,n,d) {}
	
	virtual ~PSParamgmVector4PerParticle(){};
	
};

class PSParamgmMatrix3PerParticle : public PSParamTypePerParticle<gmMatrix3>
{
public:
	virtual std::string get(int i);
	virtual void set(int i, std::string s);

	PSParamgmMatrix3PerParticle(ParticleStuff *parent, std::vector<gmMatrix3> *v,std::string sn,const std::string& n,const std::string& d = "")
	: PSParamTypePerParticle<gmMatrix3>(parent,v,sn,n,d) {}
	
	virtual ~PSParamgmMatrix3PerParticle(){};
};

class PSParamButton : public PSParam
{
public:
	class Callback 
	{
		public: 
		Callback(){};
		virtual ~Callback(){};
		virtual void onbuttonpress() = 0; 
	};

	Callback *func;		///< Function to call when button is pressed

	PSParamButton(ParticleStuff *parent, Callback *callback, const std::string& sn, const std::string& n, const std::string& d)
	: PSParam(parent,sn,n,d) { func = callback; }

    ~PSParamButton() { delete func; func = NULL; }

	virtual std::string get() { return ""; }
	virtual void set(std::string) {};
};

/*old
class PSParamComboBox : public PSParam
{
public:
	class Callback 
	{
		//realize, that the callback function receives a std::string representing the choice
		//therefore it is only natural to have the callback define the choices in the combo box
		//this is why it includes the choices vector - ELMAR
	public:
		std::vector<std::string> choices;

		Callback(){};
		virtual ~Callback(){};
		virtual void itemselect(const std::string & selection) = 0; 
	};

	Callback *func;		///< Function to call when item is selected

	PSParamComboBox  (ParticleStuff *parent, Callback *callback, 
						const std::string& shortname, 
						const std::string& name, 
						const std::string& description
						);

	~PSParamComboBox() { delete func; func = NULL; }

	virtual std::string get() { return ""; }
	virtual void set(std::string) {};
};
*/
class PSParamComboBox : public PSParam
{
public:
	class Callback 
	{
		//realize, that the callback function receives a std::string representing the choice
		//therefore it is only natural to have the callback define the choices in the combo box
		//this is why it includes the choices vector - ELMAR	

	protected:
		//these are the choices the user wants they should be set in the constructor
		std::vector<std::string> _choices;
	private:
		unsigned _selected;		
	public:
		unsigned int nbChoices() const;
		const std::string & getChoiceString(unsigned int) const; 
		const std::string & getSelectedString () const;
		unsigned int getSelected () const;
		void setSelectedItem(unsigned int i);

		//the choices vector should be initialized by your class in the constructor.
		Callback();
		virtual ~Callback(){};
		//this function sets the selected value accordingly
		//DO NOT OVERWRITE, I would like to make it private, but then the
		//interface needs to know this class.
		const void itemselectFromInterface(const std::string & selection);
		
		//this is the function you want to override:
		//it gets the selection as an unsigned int position in the choices vector and you can do whatever you want
		//if you prefer the string call _choices[i];
		virtual void itemselected(unsigned int i) = 0; 
	};

	Callback *func;		///< Function to call when item is selected

	PSParamComboBox(ParticleStuff *parent, Callback *callback, 
						const std::string& shortname, 
						const std::string& name, 
						const std::string& description
						);

	~PSParamComboBox () { delete func; func = NULL; }

	virtual std::ostream &put(std::ostream &out);
	virtual std::string get();
	virtual void set(std::string itemNumber);
	virtual PSParam &operator<<(std::string val);
};



class PSParams : public std::vector<PSParam *>
{
public:
	PSParam *findparam(const std::string &sn);

	int qlen() { return (int)size(); }

	void getq(double *q) {
		for (int i = 0; i < qlen(); i++)
			q[i] = atof(this->at(i)->get().c_str());
	}

	void setq(double *q) {
		for (int i = 0; i < qlen(); i++) {
			char buf[32];
			sprintf(buf,"%g",q[i]);
			this->at(i)->set(std::string(buf));
		}
	}

	void qname(char **qn) {
		for (int i = 0; i < qlen(); i++) {
			qn[i] = new char[this->at(i)->name.length() + 1];
			strcpy(qn[i],this->at(i)->name.c_str());
		}
	}
};

std::ostream &operator<<(std::ostream &out, const PSParams &pvec);
std::istream &operator>>(std::istream &in, PSParams &pvec);

class PSParamsPerParticle : public std::vector<PSParamPerParticle *>
{
public:
	PSParamPerParticle *findparam(const std::string &sn);
};

std::ostream &operator<<(std::ostream &out, const PSParamsPerParticle &pvec);
//std::istream &operator>>(std::istream &in, PSParamsPerParticle &pvec); // not implemented yet 

/** This is an individual record indicating that an attribute should be attached to the ParticleStuff
 * in whose attachedattributes std::vector it is held.
 */
class AttachedAttribute {
public:
	ParticleStuff *stuff;
	//ParticleAttribute **attr_var;

	std::string attr_name;	///< Name of the attribute we want
	std::string attr_name_default;		///< Default name (used to avoid redundant output)

	std::string shortname;	///< descriptor of the attachment for serialization
	std::string name;		///< descriptor of the attachment for user interface
	std::string desc;		///< long descriptor for additional help

	virtual ParticleAttribute* operator()() { return NULL; }

	virtual void attach() = 0;
	
	virtual ~AttachedAttribute(){};
};

/** This is an class that one creates in a constructor that automatically adds the attribute attachment
 * into the attachedattributes list in the ParticleStuff (passed as *parent).
 */
template<class attrclass>
class Attached : public AttachedAttribute
{
public:
	attrclass **attr_var;

	ParticleAttribute* operator()() { return *attr_var; }

	Attached(ParticleStuff *parent, attrclass **av, std::string an,const std::string& sn,const std::string& n,const std::string& d = "")
	{
		this->stuff = parent;
		attr_var = av;
		*attr_var = NULL;
		attr_name = an;
		attr_name_default = an;
		shortname = sn;
		name = n;
		desc = d;
		this->stuff->attachedattributes.push_back(this);
	}

	Attached(ParticleStuff *parent, attrclass **av)
	{
		stuff = parent;
		attr_var = av;
		*attr_var = NULL;
		std::string an = attrclass::registry_name;
		an = an.substr(an.find(':')+1);	// trim off the preface Attribute:
		attr_name_default = an;
		shortname = name = attr_name = an;
		this->stuff->attachedattributes.push_back(this);
	}

	virtual void attach() { this->stuff->attachAttribute(*attr_var,attr_name); }
	
	virtual ~Attached(){};
};

/** An instance of this class maintains the list of attributes attached to a given particle stuff
 * as the member ParticleStuff::attachedattributes.
 */
class AttachedAttributes : public std::vector<AttachedAttribute *>
{
public:
	void attach() {
		for (size_t i = 0; i < size(); i++)
		{
			at(i)->attach();
//			std::string ths = name;
			std::string n = at(i)->name;
//			int removeThisLater = 0; //tk remove this
		}
	}

	AttachedAttribute *findaa(const std::string &sn);
};

std::ostream &operator<<(std::ostream &out, const AttachedAttributes &aas);
std::ostream &operator<<(std::ostream &out, const AttachedAttribute &aa);
std::istream &operator>>(std::istream &in, AttachedAttributes &aas);

#endif

