#include "RBFBasicFunctionsAndSelector.h"
//we add the selection box to the parent surface and keep a pointer to the surface's phi functions
RBFBasicFunctionSelector::RBFBasicFunctionSelector(Surface * surface, RBFBasicFunction** f)
:SurfParamComboBox::Callback()
,_basicFunction(f)
{
	new SurfParamDouble(surface,&_cParam,0,
		"square value of the c parameter in the phi function","cParamSquare",
		"For example this parameter can be used during evaluation to smooth the surface");

	//We could replace this by an own factory... but to me it seems overkill...
	_choices.push_back(std::string("Standard pow(x,3)"));
	_choices.push_back(std::string("Standard pow(x,3) with c"));

	(*_basicFunction)= new RBFBasicFunction();
}

RBFBasicFunctionSelector::~RBFBasicFunctionSelector()
{
	delete (*_basicFunction);
	(*_basicFunction)=0;
}



void RBFBasicFunctionSelector::itemselected(unsigned int selection)
{
	if (selection==0)//standard
	{
		delete (*_basicFunction);
		(*_basicFunction)=new RBFBasicFunction();
	}
	else 
	{
		delete (*_basicFunction);
		(*_basicFunction)=new RBFBasicFunctionWithC(_cParam);		
	}
}



///////////////////////////////////////////////////
//Implementation of the simplest RBF Basic Function
///////////////////////////////////////////////////
/** Evaluate the radial basis function.
 * For 3-D, the radial basis function is ||x||^3.
 */
double RBFBasicFunction::phi(const gmVector3 & t)
{
	return (pow(t.lengthSquared(),1.5));
}

/** Evaluate the radial basis function over an interval.
 * @param x The interval over with to evaluate ||x||^3
 */
Intervald RBFBasicFunction::phi(const Box<double> & b)
{
  Box3d b1 = Box3d(b[0],b[1],b[2]);
  return b1.length().pow(3);
}

gmVector3 RBFBasicFunction::gradphi(const gmVector3 & x)
{
  double r = x.length();
  return 3.0*r*x;
}

Box3d RBFBasicFunction::gradphi(const Box<double> & b)
{
  Box3d b1 = Box3d(b[0],b[1],b[2]);
  b1 *= b.length();
  b1 *= Intervald(3.0);
  return b1;
}

gmMatrix3 RBFBasicFunction::hessphi(const gmVector3 & x)
{
  double length = x.length();

  // Cheat - added 6/5/04 by Mike Flavin
  //  - If a silhouette particle is in the same location as a control particle
  //    (which it always will be, since all particles in RBFs start at the control points),
  //    the program will crash because of divide by zero here.  So, I'm cheating by making
  //    sure length is NEVER zero.  Hopefully this will work.
  if(gmIsZero(length))
	  length = 0.000001;

  return 3.0 * 
    (gmMatrix3::identity() * length + outer(x/length,x));
}

IMatrix3d RBFBasicFunction::hessphi(const Box<double> & b)
{ 
  int i,j;
  IMatrix3d bsq;
  IMatrix3d sum;
  Box3d b1 = Box3d(b[0],b[1],b[2]);

  bsq[0][0] = bsq[1][1] = bsq[2][2] = b1.lengthSquared();
  sum = IMatrix(bsq) + IMatrix(outer(b1,b1));
  for (i = 0; i < sum.num_rows(); i++)
    for (j = 0; j < sum.num_cols(); j++)
      sum[i][j] *= (Intervald(3.0) / b1.length());

  return sum;
}






/////////////////////////////////////////////////////
//Implementation of RBF Basic Function with constant
/////////////////////////////////////////////////////

RBFBasicFunctionWithC::RBFBasicFunctionWithC(const double & parameter)
: _cParam(parameter)
{

}

double RBFBasicFunctionWithC::phi(const gmVector3 & t)
{
	return (pow(t.lengthSquared()+_cParam,1.5));
}

/** Evaluate the radial basis function over an interval.
 * @param x The interval over with to evaluate ||x||^3
 */
Intervald RBFBasicFunctionWithC::phi(const Box<double> & b)
{
	//TODO
	//this should be possible to accelerate with the right functions
	//added to box... 
	Box3d b1 = Box3d(b[0],b[1],b[2]);
	return b1.lengthSquared() + Interval<double>(_cParam,_cParam).sqrt().pow(3);
}

gmVector3 RBFBasicFunctionWithC::gradphi(const gmVector3 & x)
{
	double r = sqrt((x.lengthSquared()+_cParam));
	return 3.0*r*x;
}

Box3d RBFBasicFunctionWithC::gradphi(const Box<double> & b)
{
  Box3d b1 = Box3d(b[0],b[1],b[2]);
  b1 *= b.length()+Interval<double>(_cParam,_cParam);
  b1 *= Intervald(3.0);
  return b1;
}

gmMatrix3 RBFBasicFunctionWithC::hessphi(const gmVector3 & x)
{
  double length = sqrt(x.lengthSquared()+_cParam);

  // Cheat - added 6/5/04 by Mike Flavin
  //  - If a silhouette particle is in the same location as a control particle
  //    (which it always will be, since all particles in RBFs start at the control points),
  //    the program will crash because of divide by zero here.  So, I'm cheating by making
  //    sure length is NEVER zero.  Hopefully this will work.
  if(gmIsZero(length))
	  length = 0.000001;

  return 3.0 * 
    (gmMatrix3::identity() * length + outer(x/length,x));
}

IMatrix3d RBFBasicFunctionWithC::hessphi(const Box<double> & b)
{ 
	//TODO this can also be optimized like I did for the vector case.
	//thus interval analysis can be made faster outer((b1/scaling), b1)
  int i,j;
  IMatrix3d bsq;
  IMatrix3d sum;
  Box3d b1 = Box3d(b[0],b[1],b[2]);

  bsq[0][0] = bsq[1][1] = bsq[2][2] = b1.lengthSquared()+_cParam;
  sum = IMatrix(bsq) + IMatrix(outer(b1,b1));

  Interval<double> scaling = (Intervald(3.0)/(bsq[0][0])).sqrt();
  
  for (i = 0; i < sum.num_rows(); i++)
    for (j = 0; j < sum.num_cols(); j++)
      sum[i][j] *= scaling;

  return sum;
}
