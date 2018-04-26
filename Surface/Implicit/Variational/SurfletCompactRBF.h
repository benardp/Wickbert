/** @file SurfletCompactRBF.h
 * Expansion of CompactRBF to handle the surflet-style normal constraints of Ohtake, et al...
 * A small "piece of surface" (linear or quadratic) is fit to the region of support,
 * and blended between CSRBF centers (similar to Partition of Unity and MLS, but not quite the same)
 * @author: Scott Kircher
 * @date: November 7, 2004
 */

#ifndef SURFLETCOMPACTRBF_H
#define SURFLETCOMPACTRBF_H

#include "Surface/Implicit/Variational/CompactRBF.h"
#include <fstream>

//utility
//compute the axis-aligned bounding box (AABB) of a pointset
template<class VectorT>
void boundingBox(const std::vector<VectorT> &points,VectorT &bb_min,VectorT &bb_max)
{
	if(points.size()<1) return;
	bb_min=bb_max=points[0];
	for(unsigned int c=1;c<points.size();c++)
	{
		for(int d=0;d<3;d++)
		{
			if(bb_min[d]>points[c][d]) bb_min[d]=points[c][d];
			if(bb_max[d]<points[c][d]) bb_max[d]=points[c][d];
		}
	}
}


/** Surflet abstract base class
  */
class Surflet
{
public:
	//fit to local surface
	virtual void fit(gmVector3 point,gmVector3 normal,const RBFModelerConstraints &points,const std::vector<int> &nz_indices,const std::vector<double> &nz_weights)=0;
	//evaluate
	virtual float eval(const gmVector3 &point) const = 0;

	//analytical gradient
	virtual gmVector3 gradient(const gmVector3 &point) const = 0;

	//analytical hessian
	virtual gmMatrix3 hessian(const gmVector3 &point) const = 0;

	virtual void writeOut(std::ofstream &out) const = 0;
	virtual void readIn(std::ifstream &in) = 0;
};

/** Linear surflet
  * linear local approximation of a surface...
  * Not even the LS plane, just the plane defined by the normal and position given
  */
class LinearSurflet
{
public:
	LinearSurflet() {d=0;}
	//fit to local surface
	virtual void fit(gmVector3 point,gmVector3 normal,const RBFModelerConstraints &points,const std::vector<int> &nz_indices,const std::vector<double> &nz_weights);
	//evaluate
	virtual float eval(const gmVector3 &point) const;

	//analytical gradient
	virtual gmVector3 gradient(const gmVector3 &point) const;

	//analytical hessian
	//since the gradient is constant, this is just zero
	virtual gmMatrix3 hessian(const gmVector3 &point) const
	{
		return gmMatrix3();
	}

	virtual void writeOut(std::ofstream &out) const
	{
		out<<" ls";
		out<<" "<<n[0]<<" "<<n[1]<<" "<<n[2]<<" "<<d;
	}

	virtual void readIn(std::ifstream &in)
	{
		in>>n[0]>>n[1]>>n[2]>>d;
	}

private:
	//plane equation
	gmVector3 n;
	float d;
};

/** Quadratic surflet
  * quadratic local approximation of surface...
  * defined as height field over plane used by LinearSurflet
  */
class QuadraticSurflet
{
public:
	QuadraticSurflet() {a=0;b=0;c=0;d=0;e=0;f=0;}
	//fit to local surface
	virtual void fit(gmVector3 point,gmVector3 normal,const RBFModelerConstraints &points,const std::vector<int> &nz_indices,const std::vector<double> &nz_weights);
	//evaluate
	virtual float eval(const gmVector3 &point) const;

	//analytical gradient
	virtual gmVector3 gradient(const gmVector3 &point) const;

	//analytical hessian
	virtual gmMatrix3 hessian(const gmVector3 &point) const;

	virtual void writeOut(std::ofstream &out) const
	{
		out<<" qs";
		out<<" "<<origin[0]<<" "<<origin[1]<<" "<<origin[2]<<" ";
		out<<" "<<uhat[0]<<" "<<uhat[1]<<" "<<uhat[2]<<" ";
		out<<" "<<vhat[0]<<" "<<vhat[1]<<" "<<vhat[2]<<" ";
		out<<" "<<what[0]<<" "<<what[1]<<" "<<what[2]<<" ";
		out<<" "<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<e<<" "<<f;
	}

	virtual void readIn(std::ifstream &in)
	{
		in>>origin[0]>>origin[1]>>origin[2];
		in>>uhat[0]>>uhat[1]>>uhat[2];
		in>>vhat[0]>>vhat[1]>>vhat[2];
		in>>what[0]>>what[1]>>what[2];
		in>>a>>b>>c>>d>>e>>f;
	}

private:
	gmVector3 uhat,vhat,what; //basis vectors
	gmVector3 origin; //origin point
	float a,b,c,d,e,f;

	//global->local conversion
	inline gmVector3 toLocal(const gmVector3 &x) const {return gmVector3(dot(x-origin,uhat),dot(x-origin,vhat),dot(x-origin,what));}
	inline gmVector3 toGlobal(const gmVector3 &x) const {return x[0]*uhat+x[1]*vhat+x[2]*what+origin;}
};

/** SurfletCompactRBF class
  */
class SurfletCompactRBF : public CompactRBF
{
public:
    MAKE_NAME();

	//default constructor
    SurfletCompactRBF();
	//interpolate constructor
	SurfletCompactRBF(int cont,double support,int surflet_deg,float bloffset,const std::vector<gmVector3> &positions, 
							 const std::vector<double> &interp_vals,
							 const std::vector<gmVector3> &normals);
	//destructor
	virtual ~SurfletCompactRBF() {clearSurflets();}

	/// Solve the constraint system to get the new interpolation surface.
    virtual void updateRBF(void);

	int unsigned qlen();
	void getq(double *q);
	void _setq(double *q);
	void getqname(char **qn);

#ifndef INTERVAL_EVAL_ONLY
    virtual double proc(const gmVector3 & x);
    virtual gmVector3 grad(const gmVector3 & x);
    virtual gmMatrix3 hess(const gmVector3 & x);
#endif

	/// Function to generate a new interpolating surface from a set of constraints.
	virtual void interpolate(std::vector<gmVector3> positions, 
							 std::vector<double> interp_vals,
							 std::vector<gmVector3> normals,
							 std::valarray<bool>& flexible, 
							 bool pos_changed);

	//for use by MultiscaleCSRBF (root of hierarchy should have a blevel of 1, all others are 0)
	//this hsould not be a GUI settable parameter (and it isn't)
	void setBaseLevelOffset(double blevel) {base_level_offset=blevel;}
	double getBaseLevelOffset() {return base_level_offset;}

protected:
	//surflet degree
	//valid values:
	//  1 - Local linear approximations
	//  2 - Local quadratic approximations
	int surflet_degree;

	//surflet array (together with centers, defines all the constraints for the surface)	
	std::vector<Surflet *> surflets;
	void clearSurflets(); //deletes as well as clears vector

	//value of implicit "outside" of surface
	double base_level_offset;

	/** Function to convert the modeler constraints to true Surflet RBF constraints
	 * unlike RBF and CompactRBF, SurfletCompactRBFs do not add extra centers for
	 * normal constraints, instead they fit a surflet to the surface at that point.
	 */
	void convert_constraints(); //this includes building the KDtree

	//create a surflet of the current degree
	Surflet *newSurflet()
	{
		if(surflet_degree==1) return (Surflet*)new LinearSurflet;
		else return (Surflet*)new QuadraticSurflet;
	}
};

#endif