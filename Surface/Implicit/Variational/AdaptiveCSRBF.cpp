/** @file AdaptiveCSRBF.h
 * Extension of SurfletCompactRBF to adaptively choose CSRBF centers
 * rather than have them be at every input data point. Following Ohtake, et al...
 * "3D Scattered Data Approximation with Adaptive Compactly Supported Radial Basis Functions"
 * @author: Scott Kircher
 * @date: November 28, 2004
 */


#include "AdaptiveCSRBF.h"

#include <algorithm>
#include <fstream>
#include "tnt/tnt.h"
#include "tnt/cmat.h"
#include "tnt/vec.h"
#include "tnt/lu.h"
#include "tnt/jama_qr.h"
#include "Surface/gmTNTconvert.h"

/*************************ADAPTIVECSRBF IMPLEMENTATION*************/

REGISTER_IMPLICIT(AdaptiveCSRBF,"Variational:AdaptiveCSRBF");

/** Default empty constructor. This currently builds an RBF sphere consisting of twelve 0 centers
 * at (+/-1, +/-1, +/-1)
 */
AdaptiveCSRBF::AdaptiveCSRBF()
{
	constraints.resize(12);

	constraints[0] = RBFModelerConstraint(gmVector3(1.0,0.0,0.0), gmVector3(1.0,0.0,0.0));
	constraints[1] = RBFModelerConstraint(gmVector3(0.0,1.0,0.0), gmVector3(0.0,1.0,0.0));
	constraints[2] = RBFModelerConstraint(gmVector3(0.0,0.0,1.0), gmVector3(0.0,0.0,1.0));
	constraints[3] = RBFModelerConstraint(gmVector3(-1.0,0.0,0.0), gmVector3(-1.0,0.0,0.0));
	constraints[4] = RBFModelerConstraint(gmVector3(0.0,-1.0,0.0), gmVector3(0.0,-1.0,0.0));
	constraints[5] = RBFModelerConstraint(gmVector3(0.0,0.0,-1.0), gmVector3(0.0,0.0,-1.0));
	constraints[6] = RBFModelerConstraint(gmVector3(1.0,1.0,0.0), gmVector3(1.0,1.0,0.0).normalize());
	constraints[7] = RBFModelerConstraint(gmVector3(0.0,1.0,1.0), gmVector3(0.0,1.0,1.0).normalize());
	constraints[8] = RBFModelerConstraint(gmVector3(-1.0,0.0,1.0), gmVector3(0.0,0.0,1.0).normalize());
	constraints[9] = RBFModelerConstraint(gmVector3(-1.0,-1.0,0.0), gmVector3(-1.0,-1.0,0.0).normalize());
	constraints[10] = RBFModelerConstraint(gmVector3(0.0,-1.0,-1.0), gmVector3(0.0,-1.0,-1.0).normalize());
	constraints[11] = RBFModelerConstraint(gmVector3(1.0,0.0,-1.0), gmVector3(1.0,0.0,-1.0).normalize());

	continuity = 2;
	support_radius = 0; //not used (well, actually it IS used... but in a very weird way (see the phi functions at the very end of this file)
	surflet_degree = 2;
	param_changed = true;
	sparsity=0.5;
	overlap_threshold=1.5;
	
	min_support_size=3.0;
	max_support_size=0.0;

	base_level_offset=0.0; //not used

	updateRBF();
}

/** Returns the number of parameters
 */
unsigned int AdaptiveCSRBF::qlen()
{
	return 5;
}

/**
 * Loads parameters into q.
 * @param q An array of doubles, of size qlen, into which to copy the parameter array of the surface.
 */
void AdaptiveCSRBF::getq(double *q)
{
	q[0] = static_cast<double>(continuity);
	q[1] = static_cast<double>(surflet_degree);
	q[2] = sparsity;
	q[3] = overlap_threshold;
	q[4] = min_support_size;
}

/** Stores parameters into q.
 * @param q An array of doubles, of size qlen, containing the values to store in the parameter vector.
 */
void AdaptiveCSRBF::_setq(double *q)
{
	if(static_cast<int>(q[0])!=continuity) param_changed=true;
	if(static_cast<int>(q[1])!=surflet_degree) param_changed=true;
	if(sparsity!=q[2]) param_changed=true;
	if(overlap_threshold!=q[3]) param_changed=true;
	if(min_support_size!=q[4]) param_changed=true;

	continuity = static_cast<int>(q[0]);
	surflet_degree = static_cast<int>(q[1]);
	sparsity = q[2];
	overlap_threshold = q[3];
	min_support_size=q[4];
	
	if(param_changed)
	{
		updateRBF();
		param_changed=false;
	}
}

/**
 * Returns a list of parameter names.
 * @param qn an array of size qlen of strings listing names of parameters
 */
void AdaptiveCSRBF::getqname(char **qn)
{
	qn[0] = "Continuity";
	qn[1] = "Surflet Degree";
	qn[2] = "Sparsity";
	qn[3] = "Overlap";
	qn[4] = "Min. Support Size";
}


#ifndef INTERVAL_EVAL_ONLY
/** Evaluate the variational implicit surface.
 * The surface consists of a list of radial basis functions and a bunch of surflets 
 */
double AdaptiveCSRBF::proc(const gmVector3 & x)
{
  double fVal = 0.0;
  
  std::vector<int> nz_centers;
  support_tree.radiusQuery(x, max_support_size,nz_centers);
  
  // If the point is outside of the support radius of any point, make its value huge.
  //  This will cause ParticleValueDeath (if it is present) to destroy any floater
  //  outside of the support region, eliminating the "runaway floaters" problem.
  if(nz_centers.size() == 0)
	  return 100.0;

  double total_phi=0.0;
  for(unsigned int i = 0; i < nz_centers.size(); i++)
  {
	  //find the actual approximation point
	  int ci=approximation_indices[nz_centers[i]];
	  
	  //skip this point if its local_support_size is zero (it's not a PU center):
	  //note, this is here for safety purposes, but it should never be the case,
	  //because approximation_indices holds precicesly those points for which
	  //the local support size > 0.0.
	  if(local_support_size[ci]<=0.0) continue;

	  double this_phi=phi(x - centers[ci].c,local_support_size[ci]);
	  if(surflets[ci])
	  {
		  fVal += (surflets[ci]->eval(x) + centers[ci].w) * this_phi;
	  }else
	  {
		  fVal += (centers[ci].w) * this_phi;
	  }
	  total_phi+=this_phi;
  }
  if(total_phi==0.0) return 100.0;

  //since we're now using normalized RBFs, divide the value by the total_phi
  fVal/=total_phi;

  return fVal;//+base_level_offset;
}

gmVector3 AdaptiveCSRBF::grad(const gmVector3 & x)
{
  //the gradient is complicated by the fact that we're normalizig the CSRBFs...
  //so... I'm using finite differences for now
  return Implicit::grad(x);

}

gmMatrix3 AdaptiveCSRBF::hess(const gmVector3 & x)
{
	return Implicit::hess(x); //just use finite-differences... such a pain to differentiate 
}

#endif

/** Update the surface
  */
void AdaptiveCSRBF::updateRBF(void)
{
//generate adaptive PU (partition of unity)

  //converts constraints, builds KDtree, and fits surflets
  convert_constraints();

}

/** Convert the "conceptual" constraints into actual Surflet RBF constraint points (fit surflets)
 */      
void AdaptiveCSRBF::convert_constraints()
{
	centers.clear();
	local_support_size.clear();
	local_density_weight.clear();
	clearSurflets();

	unsigned int n = constraints.size();
	  // Insert all of the centers into the kd-tree.
	std::vector<gmVector3> points(n);
	for(size_t k = 0; k < n; k++)
	{
		points[k]=constraints[k].c;
	}
	support_tree.generate(points);

	//find bounding box of data
	gmVector3 bb_min,bb_max;
	boundingBox(points,bb_min,bb_max);
	diagonal=(bb_max-bb_min).length();

	// Convert all the constraint points.  Normal constraints get surflets fitted to them
	//  offset function value constraints.
	local_density_weight.resize(points.size(),1.0);
	local_support_size.resize(points.size(),0.0);
	surflets.resize(points.size(),0);
	max_support_size = 0; //the max_support_size is used to store the MAXIMUM support radius of all the CSRBF centers

	for(unsigned int i = 0; i < constraints.size(); i++)
	{
		centers.push_back(RBFPoint(constraints[i].c, constraints[i].h));
		centers.back().w=0.0; //initially set RBF weight to zero (i.e. we're looking just at the PU portion)
	}

	//selection of approximation centers algorithm given in the Ohtake paper
	std::vector<double> overlap(constraints.size(),0.0);
	while(true)
	{
		//find all points with overlap less than the threshold
		std::vector<int> viable_points;
		for(size_t c=0;c<overlap.size();c++)
		{
			if(overlap[c]<overlap_threshold && constraints[c].nc) viable_points.push_back(c);
		}
		//if there are no viable points, we're done!
		if(viable_points.size()<=0) break;
		//randomly permute the viable points
		std::random_shuffle(viable_points.begin(),viable_points.end());
		//choose the minimum overlap point among the first 15 viable points
		double minoverlap=overlap_threshold+1.0;
		int minpoint=0;
		for(size_t c=0;c<viable_points.size() && c<15;c++)
		{
			if(minoverlap>overlap[viable_points[c]])
			{
				minoverlap=overlap[viable_points[c]];
				minpoint=viable_points[c];
			}
		}

		//prevent this point from being chosen again
		overlap[minpoint]=overlap_threshold;

		//find the optimal support size for this point
		local_support_size[minpoint]=findOptimalSupportSize(constraints[minpoint].c,constraints[minpoint].n);
		if(local_support_size[minpoint]<min_support_size) local_support_size[minpoint]=min_support_size;

		//add a PU center here
		surflets[minpoint]=newSurflet();
		//find weights of all points near us
		std::vector<int> nz_centers;
		std::vector<double> nz_weights;
		computeLocalWeights(constraints[minpoint].c,local_support_size[minpoint],nz_centers,nz_weights);
		surflets[minpoint]->fit(constraints[minpoint].c,constraints[minpoint].n,constraints,nz_centers,nz_weights);

		//keep track of the maximum support size
		if(local_support_size[minpoint]>max_support_size) max_support_size=local_support_size[minpoint];
		
		//update the overlaps of all points
		for(size_t c=0;c<overlap.size();c++)
		{
			if(c!=minpoint)
			{
				overlap[c]+=phi(constraints[c].c-constraints[minpoint].c,local_support_size[minpoint]);
			}
		}
	}

	//rebuild the KD-tree to include ONLY the actual approximation centers.
	//if we include all the points, our radius queries later will be way slow.
	points.clear();
	approximation_indices.clear();
	for(size_t k = 0; k < n; k++)
	{
		if(local_support_size[k]>0.0)
		{
			approximation_indices.push_back(k);
			points.push_back(constraints[k].c);
		}
	}
	support_tree.generate(points);


	//output the approximation centers, so that we don't have to compute this again (its soo slow!)
	std::ofstream out("AdaptiveCSRBF_Precomp.imp");
	out<<"v0.3\n";
	out<<"Variational:AdaptiveCSRBF_Precomp "<<continuity<<" "<<approximation_indices.size();	
	for(size_t c=0;c<approximation_indices.size();c++)
	{
		int i=approximation_indices[c];
		out<<" "<<constraints[i].c[0]<<" "<<constraints[i].c[1]<<" "<<constraints[i].c[2];
		out<<" "<<local_support_size[i];
		if(surflets[i]) surflets[i]->writeOut(out);
		else out<<" ns";
	}
	out<<"\nend\n";
}	


//adapters for local support-size
double AdaptiveCSRBF::phi(const gmVector3 & x,double support)
{
	support_radius=support;
	return CompactRBF::phi(x); //base call
}

gmVector3 AdaptiveCSRBF::gradphi(const gmVector3 & x,double support)
{
	support_radius=support;
	return CompactRBF::gradphi(x); //base call
}

gmMatrix3 AdaptiveCSRBF::hessphi(const gmVector3 & x,double support)
{
	support_radius=support;
	return CompactRBF::hessphi(x); //base call
}

/**************OPTIMIZATION FUNCTIONS*************/
//compute the weighting of points (via a CS-RBF of the given size) around a given point
//returns the total sum of all weights (for optional normalization)
double AdaptiveCSRBF::computeLocalWeights(const gmVector3 &point,double sigma,std::vector<int> &nz_centers,std::vector<double> &nz_weights)
{
	//find all points near us
	support_tree.radiusQuery(point, sigma,nz_centers);

	//get weights
	double totalweight=0.0;
	nz_weights.resize(nz_centers.size(),1.0);
	for(size_t c=0;c<nz_centers.size();c++)
	{
		nz_weights[c]=local_density_weight[nz_centers[c]]*phi(constraints[nz_centers[c]].c-point,sigma);
		totalweight+=nz_weights[c];
	}
	return totalweight;
}

//In order to do the adaptive PU, we need to evaluate the error term
//for a given support size. This function returns exactly that (by first fitting
//a surflet with the desired support radius, and then evaluating the error metric
//(as given in the Ohtake paper)
//   point - surflet center
//   normal - surflet normal
//   sigma - support size (this is the parameter over which we are optimizing)
double AdaptiveCSRBF::sparseApproxError(const gmVector3 &point,const gmVector3 &normal,double sigma)
{
	//fit surflet with this support size
	Surflet *surflet = newSurflet();
	std::vector<int> nz_centers;
	std::vector<double> nz_weights;
	double totalweight=computeLocalWeights(point,sigma,nz_centers,nz_weights);
	surflet->fit(point,normal,constraints,nz_centers,nz_weights);

	//use same weighting scheme to compute Elocal^2 error
	double Elocal2=0.0;
	for(size_t c=0;c<nz_centers.size();c++)
	{
		double g=surflet->eval(constraints[nz_centers[c]].c);
		Elocal2+=(nz_weights[c]*g*g)/surflet->gradient(constraints[nz_centers[c]].c).lengthSquared();
	}
	Elocal2/=totalweight;
	//make scale-independent
	Elocal2/=(diagonal*diagonal);

	//cleanup
	delete surflet;
	
	//return objective function value
	return Elocal2 + (sparsity*sparsity*diagonal*diagonal)/(sigma*sigma);
}

//use golden-section search to minimize the objective function (equation 10 of the Ohtake paper)
double AdaptiveCSRBF::findOptimalSupportSize(const gmVector3 &point,const gmVector3 &normal)
{

	const double golden=(sqrt(5.0)-1.0)/2.0;
	double a = 0.0001; //initial endpoints
	double b = diagonal/2.0; //attempting to bracket minimum (major hack, we just pulled these numbers out of the air)
	double x1 = a+(1.0-golden)*(b-a);
	double x2 = a+golden*(b-a);
	double f1 = sparseApproxError(point,normal,x1);
	double f2 = sparseApproxError(point,normal,x2);
	while((b-a)>diagonal*0.0001)
	{
		if(f1 > f2)
		{
			a=x1;
			x1=x2;
			f1=f2;
			x2=a+golden*(b-a);
			f2= sparseApproxError(point,normal,x2);
		}else
		{
			b=x2;
			x2=x1;
			f2=f1;
			x1=a+(1.0-golden)*(b-a);
			f1= sparseApproxError(point,normal,x1);
		}
	}
	return (b+a)/2.0;
}
