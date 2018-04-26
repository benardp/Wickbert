/**
 * @file ProjNormImplicit.h
 *
 * @brief Description of ProjNormImplicit class.
 */

#ifndef PROJNORMIMPLICIT_H
#define PROJNORMIMPLICIT_H

#include "MLSSurface.h"

using namespace MLS;

/**
 * @brief MLS Surface represented nearly using the implicit presented by Adamson and Alexa.
 */
class ProjNormImplicit : public MLSSurface {

public:
    ProjNormImplicit();

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