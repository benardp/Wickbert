/** @file MultiscaleCSRBF.h
 * Uses a hierarchy of SurfletCompactRBFs to interpolate a set of points
 * in a multiscale manner.
 * @author: Scott Kircher
 * @date: November 7, 2004
 */

#include "MultiscaleCSRBF.h"

/****************************GRID CLUSTERING UTILITIES*************************/
#include <vector>
#include <map>

//cluster a set of points and normals using the specified number of subdivisions of the bounding box
//this function is efficient only if the occupied cells are sparse (this will be the case when used
//to cluster points making up a surface).
template<class VectorT>
void gridCluster(const std::vector<VectorT> &points,const std::vector<VectorT> &normals,const VectorT &bb_min,const VectorT &bb_max,int subd,std::vector<VectorT> &output_points,std::vector<VectorT> &output_normals)
{
	output_points.clear();
	output_normals.clear();
	//typedef std::pair<std::pair<int,VectorT>,VectorT> GridCell;
	std::map<int,std::pair<std::pair<int,VectorT>,VectorT> > grid;
	//typedef std::map<int,std::pair<std::pair<int,VectorT>,VectorT> >::iterator GridCellMapIter;
	//typedef GridCellMap::iterator GridCellMapIter;
	//GridCellMap grid;
	
	int dim=3;

	if(points.size()!=normals.size()) return; //eeehhh
	//run through all the points
	for(unsigned int c=0;c<points.size();c++)
	{
		//compute grid index from spatial position
		int coord;
		int index=0;
		int mult=1;
		for(int d=0;d<dim;d++)
		{
			coord=((int)((points[c][d]-bb_min[d])*((float)subd)/(bb_max[d]-bb_min[d])));
			if(coord>=subd) coord=subd-1; //clamp
			index+=coord*mult; //accumulate index
			mult*=subd;
		}
		//add this point to the grid cell
		std::pair<std::pair<int,VectorT>,VectorT> &gridcell=grid[index];
		gridcell.first.first++; //number of points in cell
		gridcell.first.second+=points[c]; //accumulated point position
		gridcell.second+=normals[c]; //accumulated normal
	}
	//generate output
	typename std::map<int,std::pair<std::pair<int,VectorT>,VectorT> >::iterator itr;
	for(itr=grid.begin();itr!=grid.end();++itr)
	{
		//average point positions
		output_points.push_back(itr->second.first.second/((float)itr->second.first.first));
		//normalize normals
		output_normals.push_back(itr->second.second.normalize());
	}
}


/****************************MULTISCALE CSRBF IMPLEMENTATION*******************/
using namespace std;

REGISTER_IMPLICIT(MultiscaleCSRBF,"Variational:MultiscaleCSRBF");

/** Default empty constructor. This currently builds an RBF sphere consisting of six 0 centers
 * at (+/-1, +/-1, +/-1)
 */
MultiscaleCSRBF::MultiscaleCSRBF()
{
	constraints.resize(6);

	constraints[0] = RBFModelerConstraint(gmVector3(1.0,0.0,0.0), gmVector3(1.0,0.0,0.0));
	constraints[1] = RBFModelerConstraint(gmVector3(0.0,1.0,0.0), gmVector3(0.0,1.0,0.0));
	constraints[2] = RBFModelerConstraint(gmVector3(0.0,0.0,1.0), gmVector3(0.0,0.0,1.0));
	constraints[3] = RBFModelerConstraint(gmVector3(-1.0,0.0,0.0), gmVector3(-1.0,0.0,0.0));
	constraints[4] = RBFModelerConstraint(gmVector3(0.0,-1.0,0.0), gmVector3(0.0,-1.0,0.0));
	constraints[5] = RBFModelerConstraint(gmVector3(0.0,0.0,-1.0), gmVector3(0.0,0.0,-1.0));

	continuity = 4;
	max_level = 12;
	support_scale = 3.0;
	surflet_degree = 2;
	param_changed = true;

	updateRBF();
	param_changed = false;
}

/** Clear all levels of the hierarchy
  */
void MultiscaleCSRBF::clearLevels()
{
	//delete all CompactRBF levels
	for(unsigned int c=0;c<levels.size();c++)
	{
		if(levels[c]) delete levels[c];
	}
	levels.clear();
}


/** Returns the number of parameters
 */
unsigned int MultiscaleCSRBF::qlen()
{
	return 4;
}

/**
 * Loads parameters into q.
 * @param q An array of doubles, of size qlen, into which to copy the parameter array of the surface.
 */
void MultiscaleCSRBF::getq(double *q)
{
	q[0] = static_cast<double>(continuity);
	q[1] = static_cast<double>(max_level);
	q[2] = support_scale;
	q[3] = static_cast<double>(surflet_degree);
}

/** Stores parameters into q.
 * @param q An array of doubles, of size qlen, containing the values to store in the parameter vector.
 */
void MultiscaleCSRBF::_setq(double *q)
{
	continuity = static_cast<int>(q[0]);
	max_level = static_cast<int>(q[1]);
	support_scale = (float) q[2];
	surflet_degree = static_cast<int>(q[3]);
	
	param_changed = true;

	updateRBF();
	param_changed = false;
}

/**
 * Returns a list of parameter names.
 * @param qn an array of size qlen of strings listing names of parameters
 */
void MultiscaleCSRBF::getqname(char **qn)
{
	qn[0] = "Continuity";
	qn[1] = "Max Level";
	qn[2] = "Support Scale";
	qn[3] = "Surflet Degree";
}

/** 
 * Loads control particles into RBF vector and calls updateRBF().
 * @param positions A vector of constraint positions.
 * @param interp_vals Is for signature compatability only. It is ignored (all points will lie on the surface).
 * @param normals A vector of the desired normals for normal constraints.
 * @param flexible A vector which serves no discernable purpose for RBFs.
 * @param pos_changed Have any of the constraints changed since we last solved for a surface?  If not, do not re-solve.
 */
void MultiscaleCSRBF::interpolate(std::vector<gmVector3> positions, std::vector<double> interp_vals, std::vector<gmVector3> normals, std::valarray<bool>& flexible, bool pos_changed)
{
	if(pos_changed || param_changed)
	{
		constraints.clear();
		constraints.resize(positions.size());

		for(unsigned int i = 0; i < positions.size(); i++)
		{
			// Normal constraint
			if(normals[i] != gmVector3(0.0, 0.0, 0.0))
				constraints[i] = RBFModelerConstraint(positions[i], normals[i], 0.0, true);
			// Non-normal constraint
			else
				constraints[i] = RBFModelerConstraint(positions[i], normals[i], 0.0, false);
		}
	
		updateRBF();
		param_changed = false;
	}
}


/**
  * Updates multiscale CSRBF with new parameters/points
  */
void MultiscaleCSRBF::updateRBF()
{
	//first clear everything
	clearLevels();
	//now, assemble a new pointset hierarchy from the set of constraints
	std::vector<gmVector3> finest_points;
	std::vector<gmVector3> finest_normals;
	for(unsigned int c=0;c<constraints.size();c++)
	{
		finest_points.push_back(constraints[c].c);
		finest_normals.push_back(constraints[c].n);
	}

	//find bounding box of data
	gmVector3 bb_min,bb_max;
	boundingBox(finest_points,bb_min,bb_max);
	double diagonal=(bb_max-bb_min).length();

	//estimate the support size based on bounding box diagonal
	double supportsize=diagonal*0.75*support_scale;

	vector<gmVector3> these_points;
	vector<gmVector3> these_normals;
	int subd=2; //initially divide grid into 8 octants
	bool baselevel=true;
	int level=1;
	do
	{
		//subdivide
		these_points.clear();
		these_normals.clear();
		gridCluster(finest_points,finest_normals,bb_min,bb_max,subd,these_points,these_normals);	

		//get values from previous level (each level is defined as an offsetting of all the ones before it)
		vector<double> these_values(these_points.size(),0.0);
		if(!baselevel)
		{
			for(unsigned int c=0;c<these_points.size();c++)
			{				
				//note, we can just use proc here, because our levels
				//vector always contains the surface approximation "so far",
				//so this is exactly what we want
///TO DO: implement the interval version
#ifndef INTERVAL_EVAL_ONLY
				these_values[c]=(-proc(these_points[c]));
#endif
			}
		}

		//create a new level
		levels.push_back(new SurfletCompactRBF(continuity,supportsize,surflet_degree,(baselevel?1.0f:0.0f),these_points,these_values,these_normals));

		baselevel=false;
		supportsize/=2.0; //supportsize decreases by a half every level
		subd<<=1; //subdivision increases by a factor of 2 at every level
		level++;

		//continue until deepest level has the same number of points as original (i.e. it is the finest level)
	}while(these_points.size()<finest_points.size() && level<=max_level);
}

#ifndef INTERVAL_EVAL_ONLY
/** Evaluate the variational implicit surface.
 * The surface consists of a hierarchy of surflet-style CS-RBFs
 */
double MultiscaleCSRBF::proc(const gmVector3 & x)
{
  double fVal = 0.0;

  //sum all the levels together
  for(unsigned int c=0;c<levels.size();c++)
  {
	  fVal+=levels[c]->proc(x);
  }
  //and that is it!
  return fVal;
}

gmVector3 MultiscaleCSRBF::grad(const gmVector3 & x)
{
  gmVector3 g(0.0,0.0,0.0);

  //sum all the levels together
  for(unsigned int c=0;c<levels.size();c++)
  {
	  g+=levels[c]->grad(x);
  }
  //and that is it!
  return g;
}

gmMatrix3 MultiscaleCSRBF::hess(const gmVector3 & x)
{
  gmMatrix3 g;

  //sum all the levels together
  for(unsigned int c=0;c<levels.size();c++)
  {
	  g+=levels[c]->hess(x);
  }
  //and that is it!
  return g;
}

#endif
