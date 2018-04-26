//@file PSClass.h
//@author Wen Su
//@date 2003-12-06
//Using spirit in boost library to parse and translate particle system language

#ifndef PSCLASS_H
#define PSCLASS_H

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;
//using namespace boost;

class PSClass
{
private:

	typedef enum {ATTRIBUTE, BEHAVIOR, SHADER} baseClassType;

	baseClassType classType;
	string className;
	string baseClassName;

	// parameters
	vector<string> parameterName;	// defaults to double
	vector<string> parameterDefaultValue;
	vector<string> parameterNamePP;
	vector<string> parameterPPType;	// std::vector<type>

	// attribute only
	string clear;
	string particleAdded;
	string particleRemoved;

	// behavior and shader
	vector<string> attachAttributeName;
	vector<string> attachAttributeType;

	// behavior only
	string applyForce;
	string applyConstraint;
	string integrate;
	string cleanup;

	// shader only
	string draw;
	string drawParticle;
	string drawShape;

	// saved tokens
	ostringstream tokens;
	// look ahead for parsing
	string nextToken;

	/// recursive parsing functions
	// ParticleAttribute, ParticleBehavior, ParticleShader name
	bool parseClass(istringstream &token);
	// Parameters or Parameters PerParticle
	bool parseParameters(istringstream &token);
	// example: jitterSize = 0.1; 
	bool parseParametersNonPP(istringstream &token);
	// example: Point oldPosition;
	bool parseParametersPP(istringstream &token);
	// attachAttributes
	bool parseAttachAttributes(istringstream &token);
	// applyForce, applyConstraint, integrate, cleanup
	bool parseBehavior(istringstream &token);
	// draw
	bool parseShader(istringstream &token);

public:
	bool parse(ifstream &in);
	bool tokenize(ifstream &in);
	bool generateHeader(ofstream &out);
	bool generateImplementation(ofstream &out);

};

#endif