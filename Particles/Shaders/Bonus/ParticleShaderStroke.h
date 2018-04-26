/**
Declaration of ParticleShaderStroke
@file ParticleShaderStroke.h
@author Wen Su
This class find the possible singular curves by connection neighboring singularities
*/

#ifndef PARTICLESHADERSTROKE_H
#define PARTICLESHADERSTROKE_H

#include <vector>

class Particles;
class ParticleShader;

class ParticleShaderStroke : public ParticleShader
{
	class Neighbor
	{
	public:
		// the first index is implicitly the index
		// the second index is stored here.
		// there are two neighbors
		unsigned int n1,n2;
	};

protected:
	double lineWidth,percentLength,curve,inkDarkness,lineTaper;
	bool autoUpdate;
	std::vector<Neighbor> chains;

public:
	MAKE_PARTICLESTUFF_NAME();

	// default constructor
	ParticleShaderStroke(Particles *ps=NULL);

	// calculate chains
	void findChains();
	void drawGradientField();

	// drawing a chain
	virtual void draw(int s);

	/// parameters
	int qlen();
	void getq(double *q);
	void setq(double *q);
	void qname(char **qn);
};

#endif