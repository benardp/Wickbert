/// Particle System Translator
/// Using spirit in boost library to parse and translate particle system language
/// auther Wen Su
/// date: 2003/12/01

#include <iostream>
#include <fstream>
#include "PSClass.h"

using namespace std;

int main(int argc, char* argv[])
{
	cout << "Particle System Translator: by Wen Su 2003/12/01"<<endl<<endl;
	// print the args
	for(int i=0;i<argc;i++)
		cout << argv[i] << " ";
	cout <<endl<<endl;

	if (argc<2)
	{
		cout << "usage:"<<endl;
		cout << "ParticleSystemTranslator filename" <<endl;
		exit(1);
	}

	// argv[0] is application name
	// argv[1] is the input file
	string psInputFileName(argv[1]);
	ifstream in(psInputFileName.c_str(), ios::in);
	if(!in)
	{
		cerr << psInputFileName << " could not be opened."<<endl;
		exit(1);
	}

	// .h file
	string psOutputFileNameH(psInputFileName.substr(0,psInputFileName.find_last_of(".")).append(".h"));
	ofstream outH(psOutputFileNameH.c_str(), ios::out);
	if(!outH)
	{
		cerr << psOutputFileNameH << " could not be opened."<<endl;
		exit(1);
	}

	// .cpp file
	string psOutputFileNameCpp(psInputFileName.substr(0,psInputFileName.find_last_of(".")).append(".cpp"));
	ofstream outCpp(psOutputFileNameCpp.c_str(), ios::out);
	if(!outCpp)
	{
		cerr << psOutputFileNameCpp << " could not be opened."<<endl;
		exit(1);
	}

	PSClass psClass;
	if (psClass.parse(in))
		cout << "Parsing " <<psInputFileName<< " succeeded." << endl;
	else
	{
		cout << "Parsing " <<psInputFileName<< " failed." << endl;
		exit(1);
	}

	// translate header
	if (psClass.generateHeader(outH))
		cout << "Output is saved into: " <<psOutputFileNameH<< endl;
	else
		cout << "generateHeader failed."<<endl;

	// translate implementation
	if (psClass.generateImplementation(outCpp))
		cout << "Output is saved into: " <<psOutputFileNameCpp<< endl;
	else
		cout << "generateImplementation failed."<<endl;

    return 0;
}