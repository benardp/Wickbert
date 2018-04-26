/**
 * @file AlexaImplicit.cpp
 *
 * @brief Implementation of AlexaImplicit class.
 */

#ifdef WB_USE_SFL

#include "AlexaImplicit.h"

// Registers the implicit the with factory.
REGISTER_IMPLICIT(AlexaImplicit,"MLS:AlexaImplicit");

/**
 * @brief Creates an AlexaImplicit.
 *
 * Creates an AlexaImplicit.  Sets the normal and energy functions we will use for
 * the MLS.
 */
AlexaImplicit::AlexaImplicit(void) :MLSSurface() {

    _normalF = MLSSurface::computeOptMLSDir;
    _energyF = MLSSurface::computeMLSEnergy;
}


/**
 * @brief Computes the value of the implicit at x.
 * 
 * Uses Adamson and Alexa's definition of an Implicit form of
 * an MLS Surface to define the value of the implicit at x.
 */
double AlexaImplicit::proc(const gmVector3 & x) {

    Vector norm;    //normal of the surface (using surflets)
    Point projP;    //projected location

    //ok - I totally hijack projP - I don't care a wit about what it's
    //value is here, I just need a variable to fill the byref call
    MLSSurface::computeAvgNormal(x, m_Kernel, norm, projP);

    //projP is the projected location of this point.
    MLSSurface::projectPoint(x, m_Kernel, _normalF, _energyF, projP);

    return dot(norm, (x - projP));
}

/**
 * Returns the gradient (normal) at x.  We use the
 * averages of the Surflet normals.
 */
gmVector3 AlexaImplicit::grad(const gmVector3 & x) {
    
    Vector norm;
    Point avg;  //return values

    MLSSurface::computeAvgNormal(x, m_Kernel, norm, avg);
    
    return norm;
}

#endif
