/**
 * @file AlexaImplicit.h
 *
 * @brief Description of AlexaImplicit class.
 */

#ifndef ALEXAIMPLICIT_H
#define ALEXAIMPLICIT_H

#include "MLSSurface.h"

using namespace MLS;

/**
 * @brief MLS Surface represented using the implicit presented by Adamson and Alexa.
 */
class AlexaImplicit : public MLSSurface {

public:
    AlexaImplicit();

    MAKE_NAME();
    
    virtual double proc(const gmVector3 & x);
    virtual gmVector3 grad(const gmVector3 & x);
	
private:
    /// Energy function type.
	typedef double (MLSSurface::*EnergyFuncType) (const Point&,Vector&,Kernel*);

    /// Normal (vector field) function type.
	typedef bool (MLSSurface::*NormalFuncType) (const Point&,Kernel*,Vector&,Point&);

    EnergyFuncType _energyF;    ///< Energy function we will use for the MLS.
    NormalFuncType _normalF;    ///< Vector field function we will use for the MLS.
	
};
#endif