/** @file Convolution.cpp 
 * @author John C. Hart
 * @date 10 Jan. 2006
 */

#include "Convolution.h"

REGISTER_IMPLICIT(Convolution,"UnaryOp:Convolution");

Convolution::Convolution(void)
{
	new SurfAttrRefParam(this,(ParticleAttribute **)&positions,"ParticlePosition","pos","positions",
		"Centers of implicit copies.");
}

