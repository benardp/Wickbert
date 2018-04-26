/*
@file ElmarsShaderTest.h
@author Elmar Eisemann
This class draws a disk for a particle object.
*/

#ifndef FINDFIRSTINTERSECTIONWITHIMPLICITS_H
#define FINDFIRSTINTERSECTIONWITHIMPLICITS_H

#include "ParticleShader.h"
class ParticleVisibility;
class ViewDependence;
class ImplicitInterrogator;
class AdaptiveRepulsionData;
class Implicit;
class ParticlePosition;

class FindFirstIntersectionWithImplicits: public ParticleShader
{

protected:
	class ImplicitOnLine
	{
	public:
		static double proc(double t, void * params);

		static gmVector3 start;
		static gmVector3 diff;
		static Implicit * imp;
	};

	static void findNegativeStepValue(	const gmVector3 & start, 
										const gmVector3 & stop, 
										const std::vector<Implicit *> implicits, 
										double stepsPerUnit, 
										gmVector3 * positionPositive, gmVector3 * positionNegative, 
										std::vector<Implicit*> * intersectedImplicits);

	static double recursiveRootSearch(	const gmVector3 & start, const gmVector3 & stop , 
										Implicit * implicit, 
										unsigned int nbSteps, double precision, unsigned int stepsPerPrecisionTest);


	ViewDependence* _view;
	ImplicitInterrogator* _impInt;
	ParticlePosition * _position;
	gmVector3 _destination;
	double _stepsPerUnit;
	int _recursive;
	double _precision;
	int _stepsPerPrecisionTest;
	bool _shootTestRay;
	bool _shootTestRaysThroughParticles;

public:
	MAKE_PARTICLESTUFF_NAME();

	FindFirstIntersectionWithImplicits(Particles *ps=NULL);
	virtual void attachAttributes();

	virtual void drawParticle(int i);
	virtual void drawPre();


	static void findIntersection(	const gmVector3 & start, 
									const gmVector3 & stop, 
									std::vector<Implicit*> implicits,
									std::pair<gmVector3, Implicit*> * result,
									double stepsPerUnit = 10, 
									unsigned int recursive = 30, 
									double precision=gmEPSILON,
									unsigned int stepsPerPrecisionTest=0);

	static bool findIntersection(	const gmVector3 & start, 
									const gmVector3 & stop, 
									Implicit* implicit,
									gmVector3 * result,
									double stepsPerUnit = 10, 
									unsigned int recursive = 30, 
									double precision=gmEPSILON,
									unsigned int stepsPerPrecisionTest=0);
};
#endif
