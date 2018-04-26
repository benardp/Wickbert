#ifndef RBFBASICFUNCTIONS_H
#define RBFBASICFUNCTIONS_H
#include "Surface/Box.h"
#include "Surface/IMatrix.h"
#include "Surface/SurfParam.h"
#include "Surface/Surface.h"

/** implementation of the RBFBasicFunction and the callback structure to 
add the selector to the interface. They are used by InterfaceRBF.
Author: Elmar Eisemann
*/

class RBFBasicFunction
{
public:
    virtual double phi(const gmVector3 & x);
    virtual gmVector3 gradphi(const gmVector3 & x);
    virtual gmMatrix3 hessphi(const gmVector3 & x);

    virtual Intervald phi(const Box<double>&);
    virtual Box3d gradphi(const Box<double>&);
    virtual IMatrix3d hessphi(const Box<double>&);
};

class RBFBasicFunctionWithC: public RBFBasicFunction
{
protected:
	const double & _cParam;
public:
	RBFBasicFunctionWithC(const double & c);
	virtual double phi(const gmVector3 & x);
    virtual gmVector3 gradphi(const gmVector3 & x);
    virtual gmMatrix3 hessphi(const gmVector3 & x);

    virtual Intervald phi(const Box<double>&);
    virtual Box3d gradphi(const Box<double>&);
    virtual IMatrix3d hessphi(const Box<double>&);

};


class RBFBasicFunctionSelector : public SurfParamComboBox::Callback
{
protected:
	double _cParam;
	RBFBasicFunction ** _basicFunction;
public:
	RBFBasicFunctionSelector(Surface * surface, RBFBasicFunction** rf);
	virtual ~RBFBasicFunctionSelector();
	virtual void itemselected(unsigned int selection);
};


#endif