/**
 * @file ProjNormImplicit.cpp
 *
 * @brief Implementation of ProjNormImplicit class.
 */
#ifdef WB_USE_SFL

#include "ProjNormImplicit.h"

// Registers the implicit the with factory.
REGISTER_IMPLICIT(ProjNormImplicit,"MLS:ProjNormImplicit");

/**
 * @brief Creates an ProjNormImplicit.
 *
 * Creates an ProjNormImplicit.  Sets the normal and energy functions we will use for
 * the MLS.
 */
ProjNormImplicit::ProjNormImplicit(void) : MLSSurface() {

    _normalF = MLSSurface::computeOptMLSDir;
    _energyF = MLSSurface::computeMLSEnergy;
}

/**
 * @brief Computes the value of the implicit at x.
 * 
 * Almost uses Adamson and Alexa's definition of an Implicit form of
 * an MLS Surface to define the value of the implicit at x.  Only
 * difference is we use the average normal from the projected
 * point.
 */
double ProjNormImplicit::proc(const gmVector3 & x) {

    Vector norm;        //normal of the surface (using surflets)
    Point projP, tempP; //projected location and a temp variable

    //projP is the projected location of this point.
    MLSSurface::projectPoint(x, m_Kernel, _normalF, _energyF, projP);

    //norm is the normal at projP
    MLSSurface::computeAvgNormal(projP, m_Kernel, norm, tempP);

    return dot(norm, x - projP);
}

/**
 * Returns the gradient (normal) at x.  We use the
 * averages of the Surflet normals.
 */
gmVector3 ProjNormImplicit::grad(const gmVector3 & x) {
    
    Vector norm;
    Point avg;  //return values

    MLSSurface::computeAvgNormal(x, m_Kernel, norm, avg);
    
    return norm;
}

#endif