/** @file SurfParam.h
 * Creates an array of parameters (usually strings) for the Surface class
 * @author John C. Hart
 * @date 29 Dec. 2004
 */

#ifndef SURFPARAM_H
#define SURFPARAM_H

//ignore deprecated warning
#pragma warning(disable : 4996)

#include <string>
#include <sstream>
#include <vector>

#include "libgm/gmVec3.h"

//itoa is not a standard C++ function. Here it is an implementation for
//compilers that don't know it -ms
#ifndef _MSC_VER
#include "itoa.h"
#endif

class Surface;
class Surfaces;
class Implicit;
class Particles;
class ParticleSystems;

class SurfParam
{
public:
	std::string shortname;
	std::string name;
	std::string desc;
	
	/** bounded? */
	bool bounded;

	SurfParam(std::string sn, std::string n, std::string d = "")
	: shortname(sn), name(n), desc(d) {}

	SurfParam() {}
	virtual ~SurfParam() {}

	virtual std::string get() = 0;
	virtual void set(std::string) = 0;
	virtual void attach(Surfaces *, ParticleSystems *) {}
};


template<class ValType>
class SurfParamType : public SurfParam
{
public:
	ValType *sv;
	
	/** Low value of the bound */
	ValType low;
	
	/** High value of the bound */
	ValType high;
	
	SurfParamType(ValType *v, ValType init,const std::string& sn,const std::string& n,const std::string& d = "")
		: SurfParam(sn,n,d) { sv = v; *sv = init; }
	
	SurfParamType(ValType *v, ValType init, ValType l, ValType h,const std::string& sn,const std::string& n,const std::string& d = "")
		: SurfParamType(sn,n,d) { sv = v; *sv = init; bounded = true; low = l; high = h; }
	
	virtual void *ref() { return sv; }
	
	virtual std::ostream &put(std::ostream &out) { return out << this->shortname << " = " << this->get(); }
	
	virtual SurfParamType &operator<<(std::string val) { set(val); return *this; }
	
	virtual ~SurfParamType(){};
};



class SurfParamBool : public SurfParamType<bool>
{
public:
	virtual std::string get() { if (*(this->sv)) return std::string("true"); else return std::string("false"); }
	virtual void set(std::string s) { *(this->sv) = (s[0] == 'T' || s[0] == 't'); }
	
	SurfParamBool(Surface *parent, bool *v, bool init,const std::string& sn,const std::string& n,const std::string& d = "");
	
	virtual ~SurfParamBool(){};
};

class SurfParamButton : public SurfParam
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
	
	SurfParamButton(Surface *parent, Callback *callback, const std::string& sn, const std::string& n, const std::string& d);
	
    ~SurfParamButton() { delete func; func = NULL; }
	
	virtual std::string get() { return "BUTTON"; }
	virtual void set(std::string) {};
};

/*implementation of a combo box to hold attributes -Elmar
*/
class SurfParamComboBox : public SurfParam
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

	SurfParamComboBox(Surface *parent, Callback *callback, 
						const std::string& shortname, 
						const std::string& name, 
						const std::string& description
						);

	~SurfParamComboBox() { delete func; func = NULL; }

	virtual std::string get();
	virtual void set(std::string itemNumber);
};


class SurfParamInt : public SurfParamType<int>
{
public:
	virtual std::string get() { char num[64]; itoa(*sv,num,10); return std::string(num); }
	virtual void set(std::string s) { *sv = atoi(s.c_str()); }

	SurfParamInt(Surface *parent, int *v, int init,const std::string& sn,const std::string& n,const std::string& d = "");
};

class SurfParamDouble : public SurfParamType<double>
{
public:
	virtual std::string get()
	{
		std::ostringstream out;
		out << *sv;
		return out.str();
	}
	virtual void set(std::string s)
	{
		std::istringstream in(s);
		in >> *sv;
	}
	SurfParamDouble(Surface *parent, double *v, double init, const std::string& sn, const std::string& n, const std::string& d = "");
	virtual ~SurfParamDouble() { }
};

class SurfParamgmVector3 : public SurfParamType<gmVector3>
{
public:
	virtual std::string get()
	{
		char buf[256];
		sprintf(buf,"(%g,%g,%g)",(*sv)[0],(*sv)[1],(*sv)[2]);
		return std::string(buf);
	}
	virtual void set(std::string s)
	{
		double x,y,z;
		sscanf(s.c_str(),"(%lf,%lf,%lf)",&x,&y,&z);
		sv->assign(x,y,z);
	}
	SurfParamgmVector3(Surface *parent, gmVector3 *v, gmVector3 init, const std::string& sn, const std::string& n, const std::string& d = "");
	virtual ~SurfParamgmVector3() { }
};

//added the SurfParamString class to replace the non-working SurfStringParam. -ms
class SurfParamString : public SurfParamType<std::string>
{
public:
	virtual std::string get() { return *(this->sv); }
	virtual void set(std::string s) { *(this->sv) = s; }
	
	SurfParamString(Surface *parent, std::string *v, std::string init,const std::string& sn,const std::string& n,const std::string& d = "");
	
	virtual ~SurfParamString(){};
};


/** This string class doesn't appear to work, at least not like the others, because instead of
    passing it a pointer to a string, it wants to maintain the string itself in its "v" member.
	Need to change this to std::string *v instead. -jch
*/

// I have tried with *v and it still does not work. - ms

class SurfStringParam : public SurfParam
{
public:
	std::string v;		///< value
	std::string def;	///< Default value

	SurfStringParam() {}
	SurfStringParam(Surface *parent, std::string *s, std::string sn, std::string n, std::string d);
	virtual ~SurfStringParam() {}
	

	virtual std::string get() { return v; }
	virtual void set(std::string s) {v = s;}
};

/** This parameter indicates a reference by name to an implicit object.
 */
class SurfImpRefParam : public SurfStringParam
{
public:
	Implicit **imp;

	SurfImpRefParam(Surface *parent, Implicit **imp_ref);

	SurfImpRefParam(Surface *parent, Implicit **imp_ref, std::string s,
							std::string sn, std::string n, std::string d);

	virtual void attach(Surfaces *surfs, ParticleSystems *psystems);
	virtual void set(std::string s) { v = s; }
	virtual ~SurfImpRefParam() {}
};

/** This parameter indicates a reference by name to a surface object.
 */
class SurfSurfRefParam : public SurfStringParam
{
public:
	Surface **surf;

	SurfSurfRefParam(Surface *parent, Surface **surf_ref);

	SurfSurfRefParam(Surface *parent, Surface **surf_ref, std::string s,
							std::string sn, std::string n, std::string d);

	virtual void attach(Surfaces *surfs, ParticleSystems *psystems);
	virtual void set(std::string s) { v = s; }
	virtual ~SurfSurfRefParam() {}
};

/** Reference by name a collection of Particles
 */
class SurfParRefParam : public SurfStringParam
{
public:
	Particles **p;

	SurfParRefParam(Surface *parent, Particles **par_ref, std::string s,
							std::string sn, std::string n, std::string d);

	virtual void attach(Surfaces *surfs, ParticleSystems *psystems);
	virtual void set(std::string s) { v = s; }

	virtual ~SurfParRefParam() {}
};

class ParticleAttribute;

/** References by name an Attribute of a Particles of a ParticleSystem.
 */
class SurfAttrRefParam : public SurfStringParam
{
public:
	ParticleAttribute **attr;

	SurfAttrRefParam(Surface *parent, ParticleAttribute **attr_ref, std::string s,
							std::string sn, std::string n, std::string d);

	virtual void attach(Surfaces *surfs, ParticleSystems *psystems);
	virtual void set(std::string s) { v = s; }
	
	virtual ~SurfAttrRefParam() {}
};

class SurfParams : public std::vector<SurfParam *>
{
public:
	void attach(Surfaces *surfs, ParticleSystems *psystems) {
		for (size_t i = 0; i < size(); i++)
			at(i)->attach(surfs,psystems);
	}

};

std::istream &operator>>(std::istream &in, SurfParams *s);
std::ostream &operator<<(std::ostream &out, SurfParams *s);
std::ostream &operator<<(std::ostream &out, SurfParam *param);

#endif
