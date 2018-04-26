/**
 * @file KolluriImplicit.h
 *
 * @brief Description of KolluriImplicit class.
 */

#ifndef KOLLURIIMPLICIT_H
#define KOLLURIIMPLICIT_H

#include "MLSSurface.h"

using namespace MLS;

/**
 * @brief MLS Surface represented using the implicit presented by Kolluri.
 */
class KolluriImplicit : public MLSSurface {

public:
    KolluriImplicit();

    MAKE_NAME();
    
    virtual double proc(const gmVector3 & x);
    virtual gmVector3 grad(const gmVector3 & x);

    virtual Box3d grad(const Box<double>& x);
    Interval<double> weightedPde1Box(const Box<double>& x, std::vector<DistancedIndex> neighbors, int arg);
    Interval<double> pde1Box(const Box<double>& x, std::vector<DistancedIndex> neighbors, int arg);
    Interval<double> normalPde1Box(const Box<double>& x, std::vector<DistancedIndex> neighbors, int arg);

private:
    /// Energy function type.
	typedef double (MLSSurface::*EnergyFuncType) (const Point&,Vector&,Kernel*);

    /// Normal (vector field) function type.
	typedef bool (MLSSurface::*NormalFuncType) (const Point&,Kernel*,Vector&,Point&);

    EnergyFuncType _energyF;    ///< Energy function we will use for the MLS.
    NormalFuncType _normalF;    ///< Vector field function we will use for the MLS.
};

#endif