/** @file AdaptiveCSRBF_Precomp.cpp
 * Same as AdaptiveCSRBF, except loaded file is expected to be a precomputed set of approximation
 * centers. Cannot be modified at runtime (no parameters)
 * @author: Scott Kircher
 * @date: December 14, 2004
 */

#include <fstream>

#include "AdaptiveCSRBF_Precomp.h"

/*************************ADAPTIVECSRBF_Precomp IMPLEMENTATION*************/

REGISTER_IMPLICIT(AdaptiveCSRBF_Precomp,"Variational:AdaptiveCSRBF_Precomp");

///Read in the (precomputed) approximation centers
bool AdaptiveCSRBF_Precomp::readImplicit(std::ifstream &file,bool verbose)
{
	//delete old surflets
	clearSurflets();
	int num;
	file>>continuity;
	file>>num;
	approximation_indices.resize(num);
	local_support_size.resize(num);
	surflets.resize(num);
	centers.resize(num);

	//read in approximation data
	max_support_size=0.0;
	for(unsigned int c=0;c<approximation_indices.size();c++)
	{
		approximation_indices[c]=c;
		file>>centers[c].c[0]>>centers[c].c[1]>>centers[c].c[2];
		centers[c].h=0.0;
		centers[c].w=0.0;
		file>>local_support_size[c];
		if(local_support_size[c]>max_support_size) max_support_size=local_support_size[c];

		std::string surflet_type;
		file>>surflet_type;
		if(surflet_type=="ls")
		{
			surflets[c]=(Surflet*)new LinearSurflet();
			surflets[c]->readIn(file);
		}else
		if(surflet_type=="qs")
		{
			surflets[c]=(Surflet*)new QuadraticSurflet();
			surflets[c]->readIn(file);
		}else
		{
			surflets[c]=0;
		}
	}

	//build KD-tree
	std::vector<gmVector3> points;
	for(unsigned int k = 0; k < centers.size(); k++)
	{
		points.push_back(centers[k].c);
	}
	support_tree.generate(points);

	return true;
}
