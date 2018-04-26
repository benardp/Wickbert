/**
 * Implementation of a class used to find the critical points of a feature contour
 * using GSL solvers.
 * @file ContourCriticalPointsGSL.cpp
 * @date 09/05/2005
 * @author Matei N. Stroila
 * @remarks
 */

#include "Surface/Implicit/Implicit.h"
#include "Surface/Implicit/Operator/RotatedImplicit.h"
#include "SingularityClassification.h"

#include "ContourCriticalPointsGSL.h"
#include "contractors.h"

INTERVAL_VECTOR fK(const INTERVAL_VECTOR& X, void* params);
INTERVAL_MATRIX JfK(const INTERVAL_VECTOR& X, void* params);
INTERVAL_VECTOR f5K(const INTERVAL_VECTOR& X, void* params);
INTERVAL_MATRIX Jf5K(const INTERVAL_VECTOR& X, void* params);

ContourCriticalPointsGSL::ContourCriticalPointsGSL(CPparams *params, double _NewtonTolerance)
{

	NewtonTolerance = _NewtonTolerance;
	WidthTolerance = 0.00001;
	SameTolerance = WidthTolerance;
	maxIter = 20;
	myF = params->F;
	pRF = new RotatedImplicit();
	pRF->setF(myF);
	params->pRF = pRF;
	params->pRinit = &Rinit;
	params->pRfinal = &Rfinal;
	singClassif = new SingularityClassification(pRF);

};

ContourCriticalPointsGSL::~ContourCriticalPointsGSL()
{
	delete pRF;
	delete singClassif;
}

int ContourCriticalPointsGSL::FindZeros4D(Box<double> &X, void *params, ContourType type, ClassifiedPointList& cpList, double* pstatus)
{
	Box<double> X1,X2; /* The two halves if we divide X */
	double width;      /* The width of the input box */
	double phi0, phi1, theta0, theta1;
	
	phi0 = ((struct CPparams *) params)->phi0; 
	theta0 = ((struct CPparams *) params)->theta0;
	phi1 = ((struct CPparams *) params)->phi1; 
	theta1 = ((struct CPparams *) params)->theta1;
	
	computeInitialFinalRotations(params);

	int numberCP = 0;
	/* Check if root is definitely NOT found in the current Box */
	if(NewtonSubdivisionBreak4D(X, params, type))
		return 0;
	
	/* See if the user wants to quit this box */
	//if(NewtonBreak(X))
	//	return;
	
	/* Use Newton's method on this box if it is small enough. */
	width = X.width();
	*pstatus = width;

	if (width < NewtonTolerance) 
    {
		if(solveContourCriticalPointGSL4D(X,params, type, cpList) == GSL_SUCCESS)
			return 1;
		else //to be removed, this else (for now, if no CP found here don't try again)
			return 0;
	}
	
	/* As a last resort, if the interval is small enough, call it a root! */
	if (width < WidthTolerance) 
    {
				
		std::cout << "added critical point in the box:" << std::endl;
		X.print();
		
		TNTPoint Xc(X.center());
		TNTPoint critPoint(5);
		for(unsigned int j = 0; j <3; j++){
			critPoint[j] = Xc[j];
		}
		double t = Xc[3];
		critPoint[3] =  phi1 * t + (1 - t) * phi0; 
		critPoint[4] =	theta1 * t + (1 - t) * theta0; 
		
		gmMatrix3 Rot = (1.0 - t) * Rinit + t * Rfinal; 
		pRF->setRotation(Rot);
		
		gmVector3 cp(critPoint[0],critPoint[1],critPoint[2]);
		
		if(!SameRoot(critPoint, cpList)) 
		{	
			cpList.push_back(ClassifiedPoint(critPoint, singClassif->getSingularityType(cp)));
			return 1;
		}
		else
		{
			std::cout << "Hmm, the subdivisionbreak method should have solved this!" << std::endl;
			return 0;
		} 
    } 
	/* Unsure, so subdivide */
	X.subdivide(X1,X2);
	numberCP += FindZeros4D(X1,params, type, cpList, pstatus);
	numberCP += FindZeros4D(X2,params, type, cpList, pstatus);
	return numberCP;
}

///Compute initial and final rotations

void  ContourCriticalPointsGSL::computeInitialFinalRotations(void *params)
{

	//get rotations
	double phi0, phi1, theta0, theta1;
	phi0 = ((struct CPparams *) params)->phi0; 
	theta0 = ((struct CPparams *) params)->theta0;
	phi1 = ((struct CPparams *) params)->phi1; 
	theta1 = ((struct CPparams *) params)->theta1;
	
	//transfrom to degree
	phi0 *= gmRADTODEG;
	theta0 *= gmRADTODEG;
	phi1 *= gmRADTODEG;
	theta1 *= gmRADTODEG;
	
	//compute initial rotation
	GLdouble modelview[16];
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
		glLoadIdentity();
		glRotated(theta0,0.0,1.0,0.0);
		glRotated(phi0,1.0,0.0,0.0);
		
		glGetDoublev(GL_MODELVIEW_MATRIX,modelview);
		gmMatrix4 _m0(modelview[0], modelview[4], modelview[8], modelview[12],
					 modelview[1], modelview[5], modelview[9], modelview[13],
					 modelview[2], modelview[6], modelview[10], modelview[14],
					 modelview[3], modelview[7], modelview[11], modelview[15]);
		for(unsigned int i = 0; i < 3; i++)
			for(unsigned int j = 0; j < 3; j++)
				Rinit[i][j] = _m0[i][j];
	glPopMatrix();

	//compute final rotation
	glPushMatrix();
		glLoadIdentity();
		glRotated(theta1,0.0,1.0,0.0);
		glRotated(phi1,1.0,0.0,0.0);
		
		glGetDoublev(GL_MODELVIEW_MATRIX,modelview);
		gmMatrix4 _m1(modelview[0], modelview[4], modelview[8], modelview[12],
					 modelview[1], modelview[5], modelview[9], modelview[13],
					 modelview[2], modelview[6], modelview[10], modelview[14],
					 modelview[3], modelview[7], modelview[11], modelview[15]);
		for(unsigned int i = 0; i < 3; i++)
			for(unsigned int j = 0; j < 3; j++)
				Rfinal[i][j] = _m1[i][j];
	glPopMatrix();

}

int ContourCriticalPointsGSL::FindZeros5D(Box<double> &X, void *params, ContourType type, ClassifiedPointList& cpList, double* pstatus)
{
	Box<double> X1,X2; /* The two halves if we divide X */
	double width;      /* The width of the input box */
	int numberCP = 0;
	/* Check if root is definitely NOT found in the current Box */
	if(NewtonSubdivisionBreak5D(X, params, type))
		return 0;
	
	/* See if the user wants to quit this box */
	//if(NewtonBreak(X))
	//	return;
	
	/* Use Newton's method on this box if it is small enough. */
	width = X.width();
	*pstatus = width;
	if (width < NewtonTolerance) 
    {
		if(solveContourCriticalPointGSL5D(X,params, type, cpList) == GSL_SUCCESS)
			return 1;
		else //to be removed, this else (for now, if no CP found here don't try again)
			return 0;
	}
	
	/* As a last resort, if the interval is small enough, call it a root! */
	if (width < WidthTolerance) 
    {
		std::cout << "added critical point in the box:" << std::endl;
		X.print();
		TNTPoint critPoint(X.center());
		std::cout << "the critical point is:" << critPoint << std::endl;
		if(!SameRoot(critPoint, cpList))
		{
			cpList.push_back(ClassifiedPoint(critPoint, UNCLASSIFIED));
			return 1;
		}
		else{
			std::cout << "Hmm, the subdivisionbreak method should have solved this!" << std::endl;
			return 0;
		} 
    } 
	
	

	/* Unsure, so subdivide */

	X.subdivide(X1,X2);
	numberCP += FindZeros5D(X1,params, type, cpList, pstatus);
	numberCP += FindZeros5D(X2,params, type, cpList, pstatus);
	return numberCP;
}

bool ContourCriticalPointsGSL::NewtonSubdivisionBreak4D(const Box<double> & X, void *params, ContourType type)
{
	bool retval = false;
	Interval<double> Fi;
		
	for (unsigned int i = 0 ; i < 4; i++){
		silhouetteContour4DCP_F(X, params, Fi,i);
		if (Fi.isNegative() || Fi.isPositive())
		{
			retval =  true;   // No root in the region - we can quit early
			break;
		}
	}	
	return retval;
}


bool ContourCriticalPointsGSL::NewtonSubdivisionBreak5D(const Box<double> & X, void *params, ContourType type)
{
	bool retval = false;
	
	Box<double> F(5);
	
	//pick the right function (depending on the contour type)
	if(type == SILHOUETTES)
		contourCP_F(X, params, F);
	
	if(type == SHADOWS)
		shadowContourCP_F(X, params, F);
	
	for (unsigned int i = 0 ; i < 5; i++){
		if (F[i].isNegative() || F[i].isPositive())
		{
			retval =  true;   // No root in the region - we can quit early
			break;
		}
	}
	
	return retval;
}


void ContourCriticalPointsGSL::print_state (size_t iter, gsl_multiroot_fdfsolver * s)
{
	std::cout << "iter =" << iter << std::endl;
	std::cout << "x =";
	for(size_t i = 0; i < s->x->size; i++)
	{
		std::cout << "/" << gsl_vector_get (s->x, i);
	}
	std::cout << std::endl;
	std::cout << "f =";
	for(size_t i = 0; i < s->f->size; i++)
	{
		std::cout << "/" << gsl_vector_get (s->f, i);
	}
	std::cout << std::endl;
}
 

int ContourCriticalPointsGSL::solveContourCriticalPointGSL5D(Box<double>& theBoxBounds, void *params, ContourType type, ClassifiedPointList& cpList)
{
	const gsl_multiroot_fdfsolver_type *T;
	gsl_multiroot_fdfsolver *s;
	
	int status;
	size_t iter = 0;
	
	const size_t n = 5;
	gsl_multiroot_function_fdf f;
	if(type == SILHOUETTES)
	{
		f.f = &contourCP_f; 
		f.df = &contourCP_df;
		f.fdf = &contourCP_fdf;
		f.n = n;
		f.params = params;
	}
	
	if(type == SHADOWS)
	{
		f.f = &shadowContourCP_f; 
		f.df = &shadowContourCP_df;
		f.fdf = &shadowContourCP_fdf;
		f.n = n;
		f.params = params;
	}
	
	double x_init[5] = {theBoxBounds[0].center(),theBoxBounds[1].center(), theBoxBounds[2].center(), theBoxBounds[3].center(), theBoxBounds[4].center()};
	gsl_vector *x = gsl_vector_alloc (n);
	
	
	for(unsigned int i = 0; i<5; i++)
		gsl_vector_set (x, i, x_init[i]);
	
//	choose the algorithm
	T = gsl_multiroot_fdfsolver_newton;
	
	s = gsl_multiroot_fdfsolver_alloc (T, n);
	gsl_multiroot_fdfsolver_set (s, &f, x);
		
	do
    {
		iter++;
		
		status = gsl_multiroot_fdfsolver_iterate (s);
		
		if (status)
			break;
		
		status = gsl_multiroot_test_residual (s->f, gmEPSILON);
    }
	while (status == GSL_CONTINUE && iter < maxIter);
		
	gsl_vector *root = gsl_multiroot_fdfsolver_root (s);
	
	double xrot = gsl_vector_get(root,3);
	double yrot = gsl_vector_get(root,4);
	double x0 = gsl_vector_get(root,0);
	double x1 = gsl_vector_get(root,1);
	double x2 = gsl_vector_get(root,2);	
 
//	if (!in(x0, theBoxBounds[0]) || !in(x1, theBoxBounds[1]) || !in(x2, theBoxBounds[2]) || !in(xrot, theBoxBounds[3]) || !in(yrot, theBoxBounds[4])) status = GSL_CONTINUE;
	 
	
	if(	status == GSL_SUCCESS) 
	{
		std::cout << "found critical point in the box:" << std::endl;
		theBoxBounds.print();
		TNTPoint critPoint(5);
		for(unsigned int j = 0; j <5; j++){
			critPoint[j] = gsl_vector_get (s->x, j);
		}
		
		if(!SameRoot(critPoint, cpList)) 
			cpList.push_back(ClassifiedPoint(critPoint, UNCLASSIFIED));
	}
	gsl_multiroot_fdfsolver_free (s);
	gsl_vector_free (x);
	return status;
}

int ContourCriticalPointsGSL::solveContourCriticalPointGSL4D(Box<double>& theBoxBounds, void *params, ContourType type, ClassifiedPointList& cpList)
{
	int status;
	unsigned int iter = 0;
	
	const size_t n = 4; //the dimension of the search space
	gsl_multiroot_function f;
	if(type == SILHOUETTES)
	{
		f.f = &silhouetteContour4DCP_f; 
		f.n = n;
		f.params = params;
	}
	
	if(type == SHADOWS)
	{
		f.f = &shadowContour4DCP_f; 
		f.n = n;
		f.params = params;
	}
	
	double x_init[4] = {theBoxBounds[0].center(),theBoxBounds[1].center(), theBoxBounds[2].center(), theBoxBounds[3].center()};
	double phi0, phi1, theta0, theta1;
	phi0 = ((struct CPparams *) params)->phi0; 
	theta0 = ((struct CPparams *) params)->theta0;
	phi1 = ((struct CPparams *) params)->phi1; 
	theta1 = ((struct CPparams *) params)->theta1;
	
	//std::cout << "[phi0,phi1]= " << phi0 * gmRADTODEG << "," << phi1 * gmRADTODEG << std::endl;
	//std::cout << "[theta0,theta1]= " << theta0 * gmRADTODEG << "," << theta1 * gmRADTODEG << std::endl;
	
	gsl_vector *x = gsl_vector_alloc (n);
	
	for(unsigned int i = 0; i<4; i++)
		gsl_vector_set (x, i, x_init[i]);
	
	//	choose the algorithm
	const gsl_multiroot_fsolver_type *T;
	gsl_multiroot_fsolver *s;
	T = gsl_multiroot_fsolver_hybrids;
	s = gsl_multiroot_fsolver_alloc (T, n);
	gsl_multiroot_fsolver_set (s, &f, x);
	
	do
    {
		iter++;
		status = gsl_multiroot_fsolver_iterate (s);
		if (status)
			break;
		status = gsl_multiroot_test_residual (s->f, 1.0e-7);
    }
	while (status == GSL_CONTINUE && iter < maxIter);

	/*
	std::cout << "status = " << gsl_strerror(status) << std::endl;
	double y0 = gsl_vector_get (s->f, 0);
	double y1 = gsl_vector_get (s->f, 1);
	double y2 = gsl_vector_get (s->f, 2);
	double y3 = gsl_vector_get (s->f, 3);
	double _width = theBoxBounds.width();
	std::cout << "the box width is :" << theBoxBounds.width() << std::endl;
	gsl_vector *root = gsl_multiroot_fsolver_root (s);
	double t = gsl_vector_get(root,3);
	double x0 = gsl_vector_get(root,0);
	double x1 = gsl_vector_get(root,1);
	double x2 = gsl_vector_get(root,2);	
	 */
	//if (!in(x0, theBoxBounds[0]) || !in(x1, theBoxBounds[1]) || !in(x2, theBoxBounds[2]) || !in(t, theBoxBounds[3])) status = GSL_CONTINUE;
	if(	status == GSL_SUCCESS) 
	{
		std::cout << "found critical point in the box:" << std::endl;
		theBoxBounds.print();
		TNTPoint critPoint(5);
		for(unsigned int j = 0; j <3; j++){
			critPoint[j] = gsl_vector_get (s->x, j);
		}
		double t = gsl_vector_get (s->x, 3);
		critPoint[3] =  phi1 * t + (1 - t) * phi0; 
		critPoint[4] =	theta1 * t + (1 - t) * theta0; 
		gmMatrix3 Rot = (1.0 - t) * Rinit + t * Rfinal; 
		pRF->setRotation(Rot);
		
		gmVector3 cp(critPoint[0],critPoint[1],critPoint[2]);
		
		if(!SameRoot(critPoint, cpList)) 
			cpList.push_back(ClassifiedPoint(critPoint, singClassif->getSingularityType(cp)));
	}
	gsl_multiroot_fsolver_free (s);
	gsl_vector_free (x);
	return status;
}

//////////////////////////
//Functions for GSL
//////////////////////////

int ContourCriticalPointsGSL::contourCP_F(const Box<double> & X, void *params, Box<double>& F)
{
	Implicit* Func  = ((struct CPparams *) params)->F;
	double radius = ((struct CPparams *) params)->radius;
	
	Box3d Pos(X[0],X[1],X[2]);
	Box3d Gradient(Func->grad(Pos));
	Box3d CameraPos( radius * sin(X[4]), (-1) * radius * sin(X[3]) * cos(X[4]), radius * cos(X[3]) * cos(X[4]));
	
	Box3d V(CameraPos - Pos);
	Box3d W(Func->hess(Pos) * V);
	Box3d G(Gradient[1]*W[2]- Gradient[2]*W[1],Gradient[2]*W[0]-
			Gradient[0]*W[2],Gradient[0]*W[1]- Gradient[1]*W[0]);
	
	F[0] = Func->proc(Pos);
	F[1] = dot(Gradient,V);
	F[2] = G[0];
	F[3] = G[1];
	F[4] = G[2];
	
	return GSL_SUCCESS;
	
}

int ContourCriticalPointsGSL::contourCP_DF(const Box<double> & X, void *params, IMatrix& Jacobian)
{
	Implicit* F  = ((struct CPparams *) params)->F;
	double radius = ((struct CPparams *) params)->radius;
	
	Box3d Pos(X[0],X[1],X[2]);
	Box3d Gradient(F->grad(Pos));
	Box3d CameraPos(radius * sin(X[4]), (-1) * radius * sin(X[3]) * cos(X[4]),
					radius * cos(X[3]) * cos(X[4]));
	Box3d V(CameraPos - Pos);
	Box3d W(F->hess(Pos) * V);	
	IMatrix Hessian(F->hess(Pos));
	
	Box3d DX3CameraPos(0.0, (-1) * radius * cos(X[3]) * cos(X[4]),
					   -1.0 * radius * sin(X[3]) * cos(X[4]));
	Box3d DX4CameraPos(radius * cos(X[4]), radius * sin(X[3]) * sin(X[4]),
					   -1.0 * radius * cos(X[3]) * sin(X[4]));
	
	for(unsigned int j = 0; j<3; j++) 
		Jacobian[0][j] = Gradient[j]; 
	
	Jacobian[0][3] = Jacobian[0][4] = 0.0;
	
	for(unsigned int j = 0; j<3; j++) 
		Jacobian[1][j] = W[j] - Gradient[j]; 
	
	Jacobian[1][3] = dot(Gradient,DX3CameraPos);	
	Jacobian[1][4] = dot(Gradient,DX4CameraPos);
	
	for(unsigned int i = 2; i<5; i++)
		for(unsigned int m = 0; m<3; m++)
			for(unsigned int j = 0; j<3; j++)
				for(unsigned int k = 0; k<3; k++)
					for(unsigned int l = 0; l<3; l++)
						Jacobian(i+1,m+1) += 
							epsilon[i-2][j][k] * Hessian(m+1,j+1) * Hessian(k+1,l+1) * V[l] + 
							epsilon[i-2][j][k] * Gradient[j] * F->D3F(Pos,m,k,l) * V[l] - epsilon[i-2][j][k] * Gradient[j] * Hessian(k+1,m+1); 
	
	for(unsigned int i = 2; i<5; i++)
		for(unsigned int j = 0; j<3; j++)
			for(unsigned int k = 0; k<3; k++)
				for(unsigned int l = 0; l<3; l++){
					Jacobian(i+1,4) += epsilon[i-2][j][k] * Gradient[j] * Hessian(k+1,l+1) * DX3CameraPos[l];
					Jacobian(i+1,5) += epsilon[i-2][j][k] * Gradient[j] * Hessian(k+1,l+1) * DX4CameraPos[l];
				}
					
					return GSL_SUCCESS;
}

int ContourCriticalPointsGSL::shadowContourCP_DF(const Box<double> & X, void *params, IMatrix& Jacobian)
{
	Implicit* F  = ((struct CPparams *) params)->F;
	double radius = ((struct CPparams *) params)->radius;
	double Lx = ((struct CPparams *) params)->Lx;
	double Ly = ((struct CPparams *) params)->Ly;
	double Lz = ((struct CPparams *) params)->Lz;
	
	Box3d Pos(X[0],X[1],X[2]);
	Box3d Gradient(F->grad(Pos));
	Box3d LightPos(
				   Lx * cos(X[4]) - Ly * sin(X[4]),
				   Ly * cos(X[3]) + Lx * sin(X[3]) * sin(X[4]) + (Lz + radius) * sin(X[3]) * cos(X[4]),
				   -1 * Ly * sin(X[3]) + Lx * cos(X[3]) * sin(X[4]) +  (Lz + radius) * cos(X[3]) * cos(X[4])
				   );
	
	Box3d DX3LightPos(
					  0.0, 
					  -1 * Ly * sin(X[3]) + Lx * cos(X[3]) * sin(X[4]) + (Lz + radius) * cos(X[3]) * cos(X[4]), 
					  -1 * Ly * cos(X[3]) - Lx * sin(X[3]) * sin(X[4]) -  (Lz + radius) * sin(X[3]) * cos(X[4])
					  );
	
	Box3d DX4LightPos(
					  -1 * Lx * sin(X[4]) - Ly * cos(X[4]),
					  Lx * sin(X[3]) * cos(X[4]) - (Lz + radius) * sin(X[3]) * sin(X[4]),
					  Lx * cos(X[3]) * cos(X[4]) -  (Lz + radius) * cos(X[3]) * sin(X[4])
					  );
	
	
	Box3d V(LightPos);
	Box3d W(F->hess(Pos) * V);	
	IMatrix Hessian(F->hess(Pos));
	
	
	for(unsigned int j = 0; j<3; j++) 
		Jacobian[0][j] = Gradient[j]; 
	
	Jacobian[0][3] = Jacobian[0][4] = 0.0;
	
	for(unsigned int j = 0; j<3; j++) 
		Jacobian[1][j] = W[j]; 
	
	Jacobian[1][3] = dot(Gradient,DX3LightPos);	
	Jacobian[1][4] = dot(Gradient,DX4LightPos);
	
	for(unsigned int i = 2; i<5; i++)
		for(unsigned int m = 0; m<3; m++)
			for(unsigned int j = 0; j<3; j++)
				for(unsigned int k = 0; k<3; k++)
					for(unsigned int l = 0; l<3; l++)
						Jacobian(i+1,m+1) += 
							epsilon[i-2][j][k] * Hessian(m+1,j+1) * Hessian(k+1,l+1) * V[l] + 
							epsilon[i-2][j][k] * Gradient[j] * F->D3F(Pos,m,k,l) * V[l]; 
	
	for(unsigned int i = 2; i<5; i++)
		for(unsigned int j = 0; j<3; j++)
			for(unsigned int k = 0; k<3; k++)
				for(unsigned int l = 0; l<3; l++){
					Jacobian(i+1,4) += epsilon[i-2][j][k] * Gradient[j] * Hessian(k+1,l+1) * DX3LightPos[l];
					Jacobian(i+1,5) += epsilon[i-2][j][k] * Gradient[j] * Hessian(k+1,l+1) * DX4LightPos[l];
				}
					
					return GSL_SUCCESS;
}




int
ContourCriticalPointsGSL::shadowContourCP_f (const gsl_vector * x, void *params, 
			 gsl_vector * f)
{
	Implicit* F  = ((struct CPparams *) params)->F;
	double radius = ((struct CPparams *) params)->radius;
	double Lx = ((struct CPparams *) params)->Lx;
	double Ly = ((struct CPparams *) params)->Ly;
	double Lz = ((struct CPparams *) params)->Lz;
	
	const double x0 = gsl_vector_get (x, 0);
	const double x1 = gsl_vector_get (x, 1);
	const double x2 = gsl_vector_get (x, 2);
	const double x3 = gsl_vector_get (x, 3);
	const double x4 = gsl_vector_get (x, 4);
	gmVector3 pos(x0,x1,x2);
	gmVector3 gradient(F->grad(pos));
	
	//transform to world coordinates
	gmVector3 lightPos(
	Lx * cos(x4) - Ly * sin(x4),
	Ly * cos(x3) + Lx * sin(x3) * sin(x4) + (Lz + radius) * sin(x3) * cos(x4),
	-1 * Ly * sin(x3) + Lx * cos(x3) * sin(x4) +  (Lz + radius) * cos(x3) * cos(x4)
					   );
	//directional light (to do: lightPos - pos for localized light
	gmVector3 v(lightPos);
	gmVector3 w(F->hess(pos) * v);
	gmVector3 g(gradient[1]*w[2]- gradient[2]*w[1],gradient[2]*w[0]-
				gradient[0]*w[2],gradient[0]*w[1]- gradient[1]*w[0]);
	
	const double y0 = F->proc(pos);
	const double y1 = dot(gradient,v);
	const double y2 = g[0];
	const double y3 = g[1];
	const double y4 = g[2];
	
	gsl_vector_set (f, 0, y0);
	gsl_vector_set (f, 1, y1);
	gsl_vector_set (f, 2, y2);
	gsl_vector_set (f, 3, y3);
	gsl_vector_set (f, 4, y4);
	
	return GSL_SUCCESS;
}

int ContourCriticalPointsGSL::shadowContourCP_F(const Box<double> & X, void *params, Box<double>& F)
{
	Implicit* Func  = ((struct CPparams *) params)->F;
	double radius = ((struct CPparams *) params)->radius;
	double Lx = ((struct CPparams *) params)->Lx;
	double Ly = ((struct CPparams *) params)->Ly;
	double Lz = ((struct CPparams *) params)->Lz;
	
	Box3d Pos(X[0],X[1],X[2]);
	Box3d Gradient(Func->grad(Pos));
	
	Box3d LightPos(
					   Lx * cos(X[4]) - Ly * sin(X[4]),
					   Ly * cos(X[3]) + Lx * sin(X[3]) * sin(X[4]) + (Lz + radius) * sin(X[3]) * cos(X[4]),
					   -1 * Ly * sin(X[3]) + Lx * cos(X[3]) * sin(X[4]) +  (Lz + radius) * cos(X[3]) * cos(X[4])
					   );
	
	Box3d V(LightPos);
	Box3d W(Func->hess(Pos) * V);
	Box3d G(Gradient[1]*W[2]- Gradient[2]*W[1],Gradient[2]*W[0]-
			Gradient[0]*W[2],Gradient[0]*W[1]- Gradient[1]*W[0]);
	
	F[0] = Func->proc(Pos);
	F[1] = dot(Gradient,V);
	F[2] = G[0];
	F[3] = G[1];
	F[4] = G[2];
	
	return GSL_SUCCESS;
	
}

int
ContourCriticalPointsGSL::shadowContourCP_df (const gsl_vector * x, void *params, 
			  gsl_matrix * J)
{
	Implicit* F  = ((struct CPparams *) params)->F;
	double radius = ((struct CPparams *) params)->radius;
	double Lx = ((struct CPparams *) params)->Lx;
	double Ly = ((struct CPparams *) params)->Ly;
	double Lz = ((struct CPparams *) params)->Lz;
	
	const double x0 = gsl_vector_get (x, 0);
	const double x1 = gsl_vector_get (x, 1);
	const double x2 = gsl_vector_get (x, 2);
	const double x3 = gsl_vector_get (x, 3);
	const double x4 = gsl_vector_get (x, 4);
	
	gmVector3 pos(x0,x1,x2);
	gmVector3 gradient(F->grad(pos));
	gmMatrix3 hessian(F->hess(pos));
	
	//transform to world coordinates
	gmVector3 lightPos(
					   Lx * cos(x4) - Ly * sin(x4),
					   Ly * cos(x3) + Lx * sin(x3) * sin(x4) + (Lz + radius) * sin(x3) * cos(x4),
					   -1 * Ly * sin(x3) + Lx * cos(x3) * sin(x4) +  (Lz + radius) * cos(x3) * cos(x4)
					   );	
	
	gmVector3 dx3LightPos(
						  0.0, 
						  -1 * Ly * sin(x3) + Lx * cos(x3) * sin(x4) + (Lz + radius) * cos(x3) * cos(x4), 
						 -1 * Ly * cos(x3) - Lx * sin(x3) * sin(x4) -  (Lz + radius) * sin(x3) * cos(x4)
						  );
	
	gmVector3 dx4LightPos(
						  -1 * Lx * sin(x4) - Ly * cos(x4),
						Lx * sin(x3) * cos(x4) - (Lz + radius) * sin(x3) * sin(x4),
						Lx * cos(x3) * cos(x4) -  (Lz + radius) * cos(x3) * sin(x4)
						  );
	
	
	gmVector3 v(lightPos);
	gmVector3 w(hessian * v);
	
	double Jacobian[5][5];
	
	for(unsigned int j = 0; j<3; j++) 
		Jacobian[0][j] = gradient[j]; 
	
	Jacobian[0][3] = Jacobian[0][4] = 0.0;
	
	for(unsigned int j = 0; j<3; j++) 
		Jacobian[1][j] = w[j]; 
	
	Jacobian[1][3] = dot(gradient,dx3LightPos);	
	Jacobian[1][4] = dot(gradient,dx4LightPos);
	
	for(unsigned int i = 2; i<5; i++)
		for(unsigned int m = 0; m<3; m++)
			for(unsigned int j = 0; j<3; j++)
				for(unsigned int k = 0; k<3; k++)
					for(unsigned int l = 0; l<3; l++)
						Jacobian[i][m] += 
							epsilon[i-2][j][k] * hessian[m][j] * hessian[k][l] * v[l] + 
							epsilon[i-2][j][k] * gradient[j] * F->D3F(pos,m,k,l) * v[l]; 
	
	for(unsigned int i = 2; i<5; i++)
		for(unsigned int j = 0; j<3; j++)
			for(unsigned int k = 0; k<3; k++)
				for(unsigned int l = 0; l<3; l++){
					Jacobian[i][3] += epsilon[i-2][j][k] * gradient[j] * hessian[k][l] * dx3LightPos[l];
					Jacobian[i][4] += epsilon[i-2][j][k] * gradient[j] * hessian[k][l] * dx4LightPos[l];
				}
					
					for(unsigned int i = 0; i<5; i++)
						for(unsigned int j = 0; j<5; j++)
							gsl_matrix_set (J, i, j, Jacobian[i][j]);
	
	return GSL_SUCCESS;
}

bool ContourCriticalPointsGSL::SameRoot(TNTPoint& newroot, ClassifiedPointList& cpList)
{
	bool same;
	ClassifiedPointList::const_iterator cpIter;
	for (cpIter = cpList.begin(); cpIter != cpList.end(); ++cpIter) 
    {
		same = true;   // Assume new point is already in the list of roots
		for (int d = 0; d < newroot.size(); d++){
			if (fabs(((*cpIter).first).operator[](d) - newroot[d]) > SameTolerance)
			{
				same = false;  // The new point and the current root are different 
				break;
			}
		}
		if (same) return true;
    }
	return false;
}

int ContourCriticalPointsGSL::shadowContour4DCP_f (const gsl_vector * x, void *params, gsl_vector * f)
{
	Implicit* F  = ((struct CPparams *) params)->F;
	double radius = ((struct CPparams *) params)->radius;
	double Lx = ((struct CPparams *) params)->Lx;
	double Ly = ((struct CPparams *) params)->Ly;
	double Lz = ((struct CPparams *) params)->Lz;
	double phi0 = ((struct CPparams *) params)->phi0; 
	double theta0 = ((struct CPparams *) params)->theta0;
	double phi1 = ((struct CPparams *) params)->phi1; 
	double theta1 = ((struct CPparams *) params)->theta1;
	
	
	const double x0 = gsl_vector_get (x, 0);
	const double x1 = gsl_vector_get (x, 1);
	const double x2 = gsl_vector_get (x, 2);
	const double t = gsl_vector_get (x, 3);
	
	double phi_t = phi1 * t + (1 - t) * phi0; //former x3
	double theta_t = theta1 * t + (1 - t) * theta0; //former x4 
	
	gmVector3 pos(x0,x1,x2);
	gmVector3 gradient(F->grad(pos));
	
	//transform to world coordinates
	gmVector3 lightPos_t(
						 Lx * cos(theta_t) - Ly * sin(theta_t),
						 Ly * cos(phi_t) + Lx * sin(phi_t) * sin(theta_t) + (Lz + radius) * sin(phi_t) * cos(theta_t),
						 -1 * Ly * sin(phi_t) + Lx * cos(phi_t) * sin(theta_t) +  (Lz + radius) * cos(phi_t) * cos(theta_t)
						 );		
	
	//directional light (to do: lightPos - pos for localized light)
	gmVector3 v(lightPos_t);
	gmVector3 w(F->hess(pos) * v);
	gmVector3 g(gradient[1]*w[2]- gradient[2]*w[1],gradient[2]*w[0]-
				gradient[0]*w[2],gradient[0]*w[1]- gradient[1]*w[0]);
	
	const double y0 = F->proc(pos) * F->proc(pos) + dot(gradient,v) * dot(gradient,v);
	const double y1 = g[0];
	const double y2 = g[1];
	const double y3 = g[2];
	
	gsl_vector_set (f, 0, y0);
	gsl_vector_set (f, 1, y1);
	gsl_vector_set (f, 2, y2);
	gsl_vector_set (f, 3, y3);

	return 0;
}



int ContourCriticalPointsGSL::shadowContour4DCP_df (const gsl_vector * x, void *params, 
					gsl_matrix * J)
{
	Implicit* F  = ((struct CPparams *) params)->F;
	double radius = ((struct CPparams *) params)->radius;
	double Lx = ((struct CPparams *) params)->Lx;
	double Ly = ((struct CPparams *) params)->Ly;
	double Lz = ((struct CPparams *) params)->Lz;
	double phi0 = ((struct CPparams *) params)->phi0; 
	double theta0 = ((struct CPparams *) params)->theta0;
	double phi1 = ((struct CPparams *) params)->phi1; 
	double theta1 = ((struct CPparams *) params)->theta1;
	
	const double x0 = gsl_vector_get (x, 0);
	const double x1 = gsl_vector_get (x, 1);
	const double x2 = gsl_vector_get (x, 2);
	const double t = gsl_vector_get (x, 3);
	
	double phi_t = phi1 * t + (1 - t) * phi0; //former x3
	double theta_t = theta1 * t + (1 - t) * theta0; //former x4 
	double phi_dt = phi1 - phi0;
	double theta_dt = theta1 - theta0;
	
	gmVector3 pos(x0,x1,x2);
	double Fpos = F->proc(pos); 
	gmVector3 gradient(F->grad(pos));
	gmMatrix3 hessian(F->hess(pos));
	
	//transform to world coordinates
	gmVector3 lightPos_t(
					   Lx * cos(theta_t) - Ly * sin(theta_t),
					   Ly * cos(phi_t) + Lx * sin(phi_t) * sin(theta_t) + (Lz + radius) * sin(phi_t) * cos(theta_t),
					   -1 * Ly * sin(phi_t) + Lx * cos(phi_t) * sin(theta_t) +  (Lz + radius) * cos(phi_t) * cos(theta_t)
					   );	
	
	gmVector3 lightPos_dt(
						  -1 * Lx * sin(theta_t) * theta_dt - Ly * cos(theta_t) * theta_dt, 
						 -1 * Ly * sin(phi_t) * phi_dt + Lx * cos(phi_t) * sin(theta_t) * phi_dt
						  + Lx * sin(phi_t) * cos(theta_t) * theta_dt
						  + (Lz + radius) * cos(phi_t) * cos(theta_t) * phi_dt
						   - (Lz + radius) * sin(phi_t) * sin(theta_t) * theta_dt, 
						  -1 * Ly * cos(phi_t) * phi_dt
						  - Lx * sin(phi_t) * sin(theta_t) * phi_dt
						  + Lx * cos(phi_t) * cos(theta_t) * theta_dt
						  -  (Lz + radius) * sin(phi_t) * cos(theta_t) * phi_dt
						  -  (Lz + radius) * cos(phi_t) * sin(theta_t) * theta_dt
						  );

	
	gmVector3 v(lightPos_t);
	double Fdotv = dot(gradient,v);
	gmVector3 w(hessian * v);
	
	double Jacobian[4][4];
	
	for(unsigned int j = 0; j<3; j++) 
		Jacobian[0][j] = 2 * Fpos * gradient[j] + 2 * Fdotv * w[j]; 
	
	Jacobian[0][3] =  2 * Fdotv * dot(gradient,lightPos_dt);
	
	
	for(unsigned int i = 1; i<4; i++)
		for(unsigned int m = 0; m<3; m++)
			for(unsigned int j = 0; j<3; j++)
				for(unsigned int k = 0; k<3; k++)
					for(unsigned int l = 0; l<3; l++)
						Jacobian[i][m] += 
							epsilon[i-1][j][k] * hessian[m][j] * hessian[k][l] * v[l] + 
							epsilon[i-1][j][k] * gradient[j] * F->D3F(pos,m,k,l) * v[l]; 
	
	for(unsigned int i = 1; i<4; i++)
		for(unsigned int j = 0; j<3; j++)
			for(unsigned int k = 0; k<3; k++)
				for(unsigned int l = 0; l<3; l++){
					Jacobian[i][3] += epsilon[i-1][j][k] * gradient[j] * hessian[k][l] * lightPos_dt[l];
				}
					
					for(unsigned int i = 0; i<4; i++)
						for(unsigned int j = 0; j<4; j++)
							gsl_matrix_set (J, i, j, Jacobian[i][j]);
	
	return 0;
}

int ContourCriticalPointsGSL::shadowContour4DCP_F(const Box<double> & X, void *params, Box<double>& F)
{
	Implicit* Func  = ((struct CPparams *) params)->F;
	double radius = ((struct CPparams *) params)->radius;
	double Lx = ((struct CPparams *) params)->Lx;
	double Ly = ((struct CPparams *) params)->Ly;
	double Lz = ((struct CPparams *) params)->Lz;
	double phi0 = ((struct CPparams *) params)->phi0; 
	double theta0 = ((struct CPparams *) params)->theta0;
	double phi1 = ((struct CPparams *) params)->phi1; 
	double theta1 = ((struct CPparams *) params)->theta1;
	
	
	Intervald T(X[3]);
	Intervald Phi_T(phi1 * T + (1.0 - T) * phi0); //former x3
	Intervald Theta_T(theta1 * T + (1.0 - T) * theta0); //former x4 

	Box3d Pos(X[0],X[1],X[2]);
	Box3d Gradient(Func->grad(Pos));
	
	//transform to world coordinates
	Box3d LightPos_T(
						 Lx * cos(Theta_T) - Ly * sin(Theta_T),
						 Ly * cos(Phi_T) + Lx * sin(Phi_T) * sin(Theta_T) + (Lz + radius) * sin(Phi_T) * cos(Theta_T),
						 -1 * Ly * sin(Phi_T) + Lx * cos(Phi_T) * sin(Theta_T) +  (Lz + radius) * cos(Phi_T) * cos(Theta_T)
						 );		
	Box3d V(LightPos_T);
	Box3d W(Func->hess(Pos) * V);
	Box3d G(Gradient[1]*W[2]- Gradient[2]*W[1],Gradient[2]*W[0]-
			Gradient[0]*W[2],Gradient[0]*W[1]- Gradient[1]*W[0]);
	
	Intervald Temp1(Func->proc(Pos));
	Intervald Temp2(dot(Gradient,V));
	F[0] = Temp1.squared() + Temp2.squared();
	F[1] = G[0];
	F[2] = G[1];
	F[3] = G[2];
	
	return 0;
	
}

int
ContourCriticalPointsGSL::contourCP_f (const gsl_vector * x, void *params, 
			 gsl_vector * f)
{
	Implicit* F  = ((struct CPparams *) params)->F;
	double radius = ((struct CPparams *) params)->radius;
	
	const double x0 = gsl_vector_get (x, 0);
	const double x1 = gsl_vector_get (x, 1);
	const double x2 = gsl_vector_get (x, 2);
	const double x3 = gsl_vector_get (x, 3);
	const double x4 = gsl_vector_get (x, 4);
	gmVector3 pos(x0,x1,x2);
	gmVector3 gradient(F->grad(pos));
	gmVector3 cameraPos(radius * sin(x4), (-1) * radius * sin(x3) * cos(x4), radius * cos(x3) * cos(x4));
	
	gmVector3 v(cameraPos - pos);
	gmVector3 w(F->hess(pos) * v);
	gmVector3 g(gradient[1]*w[2]- gradient[2]*w[1],gradient[2]*w[0]-
				gradient[0]*w[2],gradient[0]*w[1]- gradient[1]*w[0]);
	
	const double y0 = F->proc(pos);
	const double y1 = dot(gradient,v);
	const double y2 = g[0];
	const double y3 = g[1];
	const double y4 = g[2];
	
	gsl_vector_set (f, 0, y0);
	gsl_vector_set (f, 1, y1);
	gsl_vector_set (f, 2, y2);
	gsl_vector_set (f, 3, y3);
	gsl_vector_set (f, 4, y4);
	
	return GSL_SUCCESS;
}

int
ContourCriticalPointsGSL::silhouetteContour4DCP_f (const gsl_vector * x, void *params, 
			 gsl_vector * f)
{
	//get parameters
	RotatedImplicit* pRF  = ((struct CPparams *) params)->pRF;
	gmMatrix3* pRinit  = ((struct CPparams *) params)->pRinit;
	gmMatrix3* pRfinal  = ((struct CPparams *) params)->pRfinal;

	//get gsl variables
	const double x0 = gsl_vector_get (x, 0);
	const double x1 = gsl_vector_get (x, 1);
	const double x2 = gsl_vector_get (x, 2);
	const double x3 = gsl_vector_get (x, 3); //time

	gmVector3 pos(x0,x1,x2);
	gmMatrix3 R = (1.0 - x3) * (*pRinit) + x3 * (*pRfinal);
	pRF->setRotation(R);

	const double y0 = pRF->proc(pos);
	const double y1 = pRF->Fz(pos);
	const double y2 = pRF->Fzz(pos);
	const double y3 = pRF->Fx(pos) * pRF->Fyz(pos) - pRF->Fy(pos) * pRF->Fxz(pos);
	
	gsl_vector_set (f, 0, y0);
	gsl_vector_set (f, 1, y1);
	gsl_vector_set (f, 2, y2);
	gsl_vector_set (f, 3, y3);
	
	return GSL_SUCCESS;
}

void ContourCriticalPointsGSL::silhouetteContour4DCP_F(const Box<double> & X, void *params, Interval<double>& Fi, const int i)
{

	RotatedImplicit* pRF  = ((struct CPparams *) params)->pRF;
	gmMatrix3* pRinit  = ((struct CPparams *) params)->pRinit;
	gmMatrix3* pRfinal  = ((struct CPparams *) params)->pRfinal;

	Intervald T(X[3]);
	double x3 = T.center();

	gmMatrix3 R = (1.0 - x3) * (*pRinit) + x3 * (*pRfinal);
	pRF->setRotation(R);

	Box<double> Pos(3);
	for(unsigned int j = 0; j < 3; j++)
		Pos[j] = X[j];

	switch (i){
		case 0:
			Fi = pRF->proc(Pos);
			break;
		case 1:
			Fi = pRF->Fz(Pos);
		    break;
		case 2:
			Fi = pRF->Fzz(Pos);
			break;
		case 3:
			Fi = pRF->Fx(Pos) * pRF->Fyz(Pos) - pRF->Fy(Pos) * pRF->Fxz(Pos);
		    break;
		default:
			std::cout << "invalid i = " << i << std::endl;
	}
}



int
ContourCriticalPointsGSL::contourCP_df (const gsl_vector * x, void *params, 
			  gsl_matrix * J)
{
	Implicit* F  = ((struct CPparams *) params)->F;
	double radius = ((struct CPparams *) params)->radius;
	
	const double x0 = gsl_vector_get (x, 0);
	const double x1 = gsl_vector_get (x, 1);
	const double x2 = gsl_vector_get (x, 2);
	const double x3 = gsl_vector_get (x, 3);
	const double x4 = gsl_vector_get (x, 4);
	
	gmVector3 pos(x0,x1,x2);
	gmVector3 gradient(F->grad(pos));
	gmMatrix3 hessian(F->hess(pos));
	
	gmVector3 cameraPos(radius * sin(x4), (-1) * radius * sin(x3) * cos(x4),radius * cos(x3) * cos(x4));
	
	gmVector3 dx3CameraPos(0.0, (-1) * radius * cos(x3) * cos(x4), (-1) * radius * sin(x3) * cos(x4));
	
	gmVector3 dx4CameraPos(radius * cos(x4), radius * sin(x3) * sin(x4), (-1) * radius * cos(x3) * sin(x4));
	
	
	gmVector3 v(cameraPos - pos);
	gmVector3 w(hessian * v);
	
	double Jacobian[5][5];
	
	for(unsigned int j = 0; j<3; j++) 
		Jacobian[0][j] = gradient[j]; 
	
	Jacobian[0][3] = Jacobian[0][4] = 0.0;
	
	for(unsigned int j = 0; j<3; j++) 
		Jacobian[1][j] = w[j] - gradient[j]; 
	
	Jacobian[1][3] = dot(gradient,dx3CameraPos);	
	Jacobian[1][4] = dot(gradient,dx4CameraPos);
	
	for(unsigned int i = 2; i<5; i++)
		for(unsigned int m = 0; m<3; m++)
			for(unsigned int j = 0; j<3; j++)
				for(unsigned int k = 0; k<3; k++)
					for(unsigned int l = 0; l<3; l++)
						Jacobian[i][m] += 
							epsilon[i-2][j][k] * hessian[m][j] * hessian[k][l] * v[l] + 
							epsilon[i-2][j][k] * gradient[j] * F->D3F(pos,m,k,l) * v[l] - epsilon[i-2][j][k] * gradient[j] * hessian[k][m]; 
	
	for(unsigned int i = 2; i<5; i++)
		for(unsigned int j = 0; j<3; j++)
			for(unsigned int k = 0; k<3; k++)
				for(unsigned int l = 0; l<3; l++){
					Jacobian[i][3] += epsilon[i-2][j][k] * gradient[j] * hessian[k][l] * dx3CameraPos[l];
					Jacobian[i][4] += epsilon[i-2][j][k] * gradient[j] * hessian[k][l] * dx4CameraPos[l];
				}
					
					for(unsigned int i = 0; i<5; i++)
						for(unsigned int j = 0; j<5; j++)
							gsl_matrix_set (J, i, j, Jacobian[i][j]);
	
	return GSL_SUCCESS;
}

int
ContourCriticalPointsGSL::contourCP_fdf (const gsl_vector * x, void *params,
			   gsl_vector * f, gsl_matrix * J)
{
	contourCP_f (x, params, f);
	contourCP_df (x, params, J);
	
	return GSL_SUCCESS;
}

int ContourCriticalPointsGSL::shadowContourCP_fdf (const gsl_vector * x, void *params, gsl_vector * f, gsl_matrix * J)
{
	shadowContourCP_f (x, params, f);
	shadowContourCP_df (x, params, J);
	
	return GSL_SUCCESS;
}

int ContourCriticalPointsGSL::shadowContour4DCP_fdf (const gsl_vector * x, void *params, gsl_vector * f, gsl_matrix * J)
{
	shadowContour4DCP_f (x, params, f);
	shadowContour4DCP_df (x, params, J);
	
	return 0;
}

//added for the new system:  F = Fv = kG = kR = 0	
int ContourCriticalPointsGSL::silh_kGkR_4DCP_f (const gsl_vector *x, void *params, gsl_vector *f)
{
	//get parameters
	Implicit* F  = ((struct CPparams *) params)->F;
	double radius = ((struct CPparams *) params)->radius;
	double phi0 = ((struct CPparams *) params)->phi0; 
	double theta0 = ((struct CPparams *) params)->theta0;
	double phi1 = ((struct CPparams *) params)->phi1; 
	double theta1 = ((struct CPparams *) params)->theta1;
	
	//get gsl variables
	const double x0 = gsl_vector_get (x, 0);
	const double x1 = gsl_vector_get (x, 1);
	const double x2 = gsl_vector_get (x, 2);
	const double x3 = gsl_vector_get (x, 3); //time
	
	gmVector3 pos(x0,x1,x2);
	gmVector3 gradient(F->grad(pos));
	double theta = (1 - x3) * theta0 + x3 * theta1;
	double phi = (1 - x3) * phi0 + x3 * phi1;
	
	gmVector3 cameraPos( radius * sin(theta), -1 * radius * sin(phi) * cos(theta), -1 * radius * cos(phi) * cos(theta));
	
	gmVector3 v(cameraPos - pos);
	
	const double y0 = F->proc(pos);
	const double y1 = dot(gradient,v);
	const double y2 = F->gaussianCurvature(pos); 
	const double y3 = dot(F->hess(pos)*v,v);
	
	gsl_vector_set (f, 0, y0);
	gsl_vector_set (f, 1, y1);
	gsl_vector_set (f, 2, y2);
	gsl_vector_set (f, 3, y3);
	
	return 0;
}

void ContourCriticalPointsGSL::silh_kGkR_4DCP_F(const Box<double>& X, void *params, Interval<double>& Fi, const int i)
{
	Implicit* Func  = ((struct CPparams *) params)->F;
	double radius = ((struct CPparams *) params)->radius;
	double phi0 = ((struct CPparams *) params)->phi0; 
	double theta0 = ((struct CPparams *) params)->theta0;
	double phi1 = ((struct CPparams *) params)->phi1; 
	double theta1 = ((struct CPparams *) params)->theta1;	
	
	Box3d Pos(X[0],X[1],X[2]);
	Box3d Gradient, CameraPos, V, P;
	Intervald theta, phi;

	switch (i){
		case 0:
			Fi = Func->proc(Pos);
			break;
		case 1:
			Gradient= Func->grad(Pos);
			theta = (1 - X[3]) * theta0 + X[3] * theta1;
			phi = (1 - X[3]) * phi0 + X[3] * phi1;
			CameraPos[0] = radius * sin(theta);
			CameraPos[1] = (-1) * radius * sin(phi) * cos(theta);
			CameraPos[2] = -1 * radius * cos(phi) * cos(theta);
			V = CameraPos - Pos;
			Fi = dot(Gradient,V);
		    break;
		case 2:
			Fi = Func->gaussianCurvature(Pos);
			break;
		case 3:
			theta = (1 - X[3]) * theta0 + X[3] * theta1;
			phi = (1 - X[3]) * phi0 + X[3] * phi1;
			CameraPos[0] = radius * sin(theta);
			CameraPos[1] = (-1) * radius * sin(phi) * cos(theta);
			CameraPos[2] = -1 * radius * cos(phi) * cos(theta);
			V = CameraPos - Pos;
			P = Func->hess(Pos)*V;
			Fi = dot(P,V);
		    break;
		default:
			std::cout << "invalid i = " << i << std::endl;
	}
}

int ContourCriticalPointsGSL::shadow_kGkR_4DCP_f (const gsl_vector *x, void *params, gsl_vector *f)
{
	//get parameters
	Implicit* F  = ((struct CPparams *) params)->F;
	double radius = ((struct CPparams *) params)->radius;
	double Lx = ((struct CPparams *) params)->Lx;
	double Ly = ((struct CPparams *) params)->Ly;
	double Lz = ((struct CPparams *) params)->Lz;
	double phi0 = ((struct CPparams *) params)->phi0; 
	double theta0 = ((struct CPparams *) params)->theta0;
	double phi1 = ((struct CPparams *) params)->phi1; 
	double theta1 = ((struct CPparams *) params)->theta1;
	
	//get gsl variables
	const double x0 = gsl_vector_get (x, 0);
	const double x1 = gsl_vector_get (x, 1);
	const double x2 = gsl_vector_get (x, 2);
	const double x3 = gsl_vector_get (x, 3); //time
	
	gmVector3 pos(x0,x1,x2);
	gmVector3 gradient(F->grad(pos));
	double theta = (1 - x3) * theta0 + x3 * theta1;
	double phi = (1 - x3) * phi0 + x3 * phi1;
	//transform to world coordinates
	gmVector3 lightPos(
					   Lx * cos(theta) - Lz * sin(theta),
					   Ly * cos(phi) + Lz * sin(phi) * cos(theta) + Lx * sin(phi) * sin(theta), -1 * Ly * sin(phi) + Lx * cos(phi) * sin(theta) +  Lz  * cos(phi) * cos(theta));
	//directional light (to do: lightPos - pos for localized light
	gmVector3 v(lightPos);
	
	const double y0 = F->proc(pos);
	const double y1 = dot(gradient,v);
	const double y2 = F->gaussianCurvature(pos); 
	const double y3 = dot(F->hess(pos)*v,v);
	
	gsl_vector_set (f, 0, y0);
	gsl_vector_set (f, 1, y1);
	gsl_vector_set (f, 2, y2);
	gsl_vector_set (f, 3, y3);
	
	return 0;
}

int ContourCriticalPointsGSL::shadow_kGkR_4DCP_F(const Box<double>& X, void *params, Box<double>& F)
{
	Implicit* Func  = ((struct CPparams *) params)->F;
	double radius = ((struct CPparams *) params)->radius;
	double Lx = ((struct CPparams *) params)->Lx;
	double Ly = ((struct CPparams *) params)->Ly;
	double Lz = ((struct CPparams *) params)->Lz;
	double phi0 = ((struct CPparams *) params)->phi0; 
	double theta0 = ((struct CPparams *) params)->theta0;
	double phi1 = ((struct CPparams *) params)->phi1; 
	double theta1 = ((struct CPparams *) params)->theta1;	
	
	Box3d Pos(X[0],X[1],X[2]);
	Box3d Gradient(Func->grad(Pos));
	Intervald theta = (1 - X[3]) * theta0 + X[3] * theta1;
	Intervald phi = (1 - X[3]) * phi0 + X[3] * phi1;
	Box3d LightPos(
				   Lx * cos(theta) - Lz * sin(theta),
				   Ly * cos(phi) + Lz * sin(phi) * cos(theta) + Lx * sin(phi) * sin(theta), -1 * Ly * sin(phi) + Lx * cos(phi) * sin(theta) +  Lz  * cos(phi) * cos(theta));
	
	Box3d V(LightPos);
	Box3d P = Func->hess(Pos)*V;
	Intervald Kr = dot(P,V);
	
	F[0] = Func->proc(Pos);
	F[1] = dot(Gradient,V);
	F[2] = Func->gaussianCurvature(Pos); 
	F[3] = Kr;
	
	return 0;
}


int ContourCriticalPointsGSL::FindZeros_kGkR_4D(Box<double> &X, void *params, ContourType type, ClassifiedPointList& cpList, double* pstatus)
{
	Box<double> X1,X2; /* The two halves if we divide X */
	double width;      /* The width of the input box */
			
	int numberCP = 0;
	/* Check if root is definitely NOT found in the current Box */
	if(NewtonSubdivisionBreak_kGkR_4D(X, params, type))
		return 0;
	
	/* See if the user wants to quit this box */
	//if(NewtonBreak(X))
	//	return;
	
	/* Use Newton's method on this box if it is small enough. */
	//width = X.width();
	width = customWidth(X);
	*pstatus = width;
	
	if (width < NewtonTolerance) 
    {
		if(solveContourCriticalPointGSL_kGkR_4D(X,params, type, cpList) == GSL_SUCCESS)
			return 1;
		else //to be removed, this else (for now, if no CP found here don't try again)
			return 0;
	}
	
	/* As a last resort, if the interval is small enough, call it a root! */
	if (width < WidthTolerance) 
    {
		
		std::cout << "added critical point in the box:" << std::endl;
		X.print();
		
		double phi0, phi1, theta0, theta1;
		phi0 = ((struct CPparams *) params)->phi0; 
		theta0 = ((struct CPparams *) params)->theta0;
		phi1 = ((struct CPparams *) params)->phi1; 
		theta1 = ((struct CPparams *) params)->theta1;
		
		TNTPoint Xc(X.center());
		TNTPoint critPoint(5);
		for(unsigned int j = 0; j <3; j++){
			critPoint[j] = Xc[j];
		}
		double t = Xc[3];
		critPoint[3] =  phi1 * t + (1 - t) * phi0; 
		critPoint[4] =	theta1 * t + (1 - t) * theta0; 
				
		gmVector3 cp(critPoint[0],critPoint[1],critPoint[2]);
		
		if(!SameRoot(critPoint, cpList)) 
		{	
			cpList.push_back(ClassifiedPoint(critPoint, singClassif->getSingularityType(cp)));
			return 1;
		}
		else
		{
			std::cout << "Hmm, the subdivisionbreak method should have solved this!" << std::endl;
			return 0;
		} 
    } 
	/* Unsure, so subdivide */
	
	double maxwidth = -1.0;
	int dim, idx;
	
	for (idx = 0; idx < 3; idx++) 
		if (X[idx].width() > maxwidth) 
		{
			dim = idx;
			maxwidth = X[idx].width();
		}
	
	
	double middle = (1.0 - 0.5) * X[dim].low() + 0.5 * X[dim].high();
	
	X1 = X2 = X;
	
	X1[dim] = Intervald(X[dim].low(),middle);
	X2[dim] = Intervald(middle,X[dim].high());
	
	//X.subdivide(X1,X2);
	numberCP += FindZeros_kGkR_4D(X1,params, type, cpList, pstatus);
	numberCP += FindZeros_kGkR_4D(X2,params, type, cpList, pstatus);
	return numberCP;
} 
	
int ContourCriticalPointsGSL::FindZerosIA(Box<double> &X, void *params, ContourType type, ClassifiedPointList& cpList, double* pstatus)
{
	Box<double> X1,X2; /* The two halves if we divide X */
	double width;      /* The width of the input box */
				
	int numberCP = 0;

	//X =  KrawczykContractor( fK,  JfK,  X, params);
	Box<double> Xo(X.center());
	Box<double> p(X - Xo);

	p = PrecGaussSeidelIteration(JfK(X,params), -fK(Xo,params), p);

	Box<double> temp(p+Xo);
	X.intersection(temp);
	
	width = customWidth(X);

	*pstatus = width;
	
	if (width < WidthTolerance) 
    {
		
		double phi0, phi1, theta0, theta1;
		phi0 = ((struct CPparams *) params)->phi0; 
		theta0 = ((struct CPparams *) params)->theta0;
		phi1 = ((struct CPparams *) params)->phi1; 
		theta1 = ((struct CPparams *) params)->theta1;
		
		TNTPoint Xc(X.center());
		TNTPoint critPoint(5);
		for(unsigned int j = 0; j <3; j++){
			critPoint[j] = Xc[j];
		}
		double t = Xc[3];
		critPoint[3] =  phi1 * t + (1 - t) * phi0; 
		critPoint[4] =	theta1 * t + (1 - t) * theta0; 
				
		gmVector3 cp(critPoint[0],critPoint[1],critPoint[2]);
		
		if(!SameRoot(critPoint, cpList)) 
		{	
			std::cout << "added critical point in the box:" << std::endl;
			X.print();			
			cpList.push_back(ClassifiedPoint(critPoint, singClassif->getSingularityType(cp)));
			return 1;
		}
		else
		{
			std::cout << "Same CP!" << std::endl;
			return 0;
		} 
    } 
	/* Unsure, so subdivide */
	/*
	double maxwidth = -1.0;
	int dim, idx;
	
	for (idx = 0; idx < 4; idx++) 
		if (X[idx].width() > maxwidth) 
		{
			dim = idx;
			maxwidth = X[idx].width();
		}
	
	
	double middle = (1.0 - 0.5) * X[dim].low() + 0.5 * X[dim].high();
	
	X1 = X2 = X;
	
	X1[dim] = Intervald(X[dim].low(),middle);
	X2[dim] = Intervald(middle,X[dim].high());
	*/
	X.subdivide(X1,X2);
	if(cpList.size()  > 100) return numberCP;
	numberCP += FindZerosIA(X1,params, type, cpList, pstatus);
	numberCP += FindZerosIA(X2,params, type, cpList, pstatus);
	return numberCP;
} 

int ContourCriticalPointsGSL::FindZeros5IA(Box<double> &X, void *params, ContourType type, ClassifiedPointList& cpList, double* pstatus)
{
	Box<double> X1,X2; /* The two halves if we divide X */
	double width;      /* The width of the input box */
				
	int numberCP = 0;

//	X =  KrawczykContractor( f5K,  Jf5K,  X, params);

	width = customWidth(X);
	*pstatus = width;
	
	if (width < WidthTolerance) 
	{
		TNTPoint critPoint(X.center());			
		gmVector3 cp(critPoint[0],critPoint[1],critPoint[2]);
		
		if(!SameRoot(critPoint, cpList)) 
		{	
			std::cout << "added critical point in the box:" << std::endl;
			X.print();			
			cpList.push_back(ClassifiedPoint(critPoint, singClassif->getSingularityType(cp)));
			return 1;
		}
		else
		{
			//std::cout << "Same CP!" << std::endl;
			return 0;
		} 
	} 
	/* Unsure, so subdivide */
	
	double maxwidth = -1.0;
	int dim, idx;
	
	for (idx = 0; idx < 5; idx++) 
		if (X[idx].width() > maxwidth) 
	{
		dim = idx;
		maxwidth = X[idx].width();
	}
	
	
	double middle = (1.0 - 0.5) * X[dim].low() + 0.5 * X[dim].high();
	
	X1 = X2 = X;
	
	X1[dim] = Intervald(X[dim].low(),middle);
	X2[dim] = Intervald(middle,X[dim].high());
	
	//X.subdivide(X1,X2);
	if(cpList.size()  > 100) return numberCP;
	numberCP += FindZeros5IA(X1,params, type, cpList, pstatus);
	numberCP += FindZeros5IA(X2,params, type, cpList, pstatus);
	return numberCP;
} 



bool ContourCriticalPointsGSL::NewtonSubdivisionBreak_kGkR_4D(const Box<double> & X, void *params, ContourType type)
{
	bool retval = false;
	
	Interval<double> Fi;
		
	for (unsigned int i = 0 ; i < 4; i++){
		silh_kGkR_4DCP_F(X, params, Fi,i);
		if (Fi.isNegative() || Fi.isPositive())
		{
			retval =  true;   // No root in the region - we can quit early
			break;
		}
	}	
	return retval;
}

int ContourCriticalPointsGSL::solveContourCriticalPointGSL_kGkR_4D(Box<double>& theBoxBounds, void *params, ContourType type, ClassifiedPointList& cpList)
{

	int status;
	unsigned int iter = 0;
	const size_t n = 4; //the dimension of the search space
	
#if 0
	gsl_multiroot_function f;
	if(type == SILHOUETTES)
	{
		f.f = &silh_kGkR_4DCP_f; 
		f.n = n;
		f.params = params;
	}
	
	
	if(type == SHADOWS)
	{
		f.f = &shadow_kGkR_4DCP_F; 
		f.n = n;
		f.params = params;
	}
	
	//	choose the algorithm
	const gsl_multiroot_fsolver_type *T;
	gsl_multiroot_fsolver *s;
	
#endif	
	
	
#if 1	
	gsl_multiroot_function_fdf f;
	if(type == SILHOUETTES)
	{
		f.f = &silh_kGkR_4DCP_f; 
		f.df = &silh_kGkR_4DCP_df;
		f.fdf = &silh_kGkR_4DCP_fdf;
		f.n = n;
		f.params = params;
	}
	
	const gsl_multiroot_fdfsolver_type *T;
	gsl_multiroot_fdfsolver *s;
	
#endif
	
	double x_init[4] = {theBoxBounds[0].center(),theBoxBounds[1].center(), theBoxBounds[2].center(), theBoxBounds[3].center()};
	double phi0, phi1, theta0, theta1;
	phi0 = ((struct CPparams *) params)->phi0; 
	theta0 = ((struct CPparams *) params)->theta0;
	phi1 = ((struct CPparams *) params)->phi1; 
	theta1 = ((struct CPparams *) params)->theta1;
	
	gsl_vector *x = gsl_vector_alloc (n);
	
	for(unsigned int i = 0; i<4; i++)
		gsl_vector_set (x, i, x_init[i]);
	
	//	choose the algorithm
#if 0
	T = gsl_multiroot_fsolver_dnewton;
	s = gsl_multiroot_fsolver_alloc (T, n);
	gsl_multiroot_fsolver_set (s, &f, x);
	
	do
    {
		iter++;
		status = gsl_multiroot_fsolver_iterate (s);
		if (status)
			break;
		status = gsl_multiroot_test_residual (s->f, 1.0e-7);
    }
	while (status == GSL_CONTINUE && iter < maxIter);
	gsl_vector *root = gsl_multiroot_fsolver_root (s);
#endif
	
	
	//	choose the algorithm
#if 1
	T = gsl_multiroot_fdfsolver_newton;
	s = gsl_multiroot_fdfsolver_alloc (T, n);
	gsl_multiroot_fdfsolver_set (s, &f, x);
	
	do
    {
		iter++;
		status = gsl_multiroot_fdfsolver_iterate (s);
		if (status)
			break;
		status = gsl_multiroot_test_residual (s->f, 1.0e-7);
    }
	while (status == GSL_CONTINUE && iter < maxIter);
	gsl_vector *root = gsl_multiroot_fdfsolver_root (s);
#endif
	
	
//	 std::cout << "status = " << gsl_strerror(status) << std::endl;
//	 double y0 = gsl_vector_get (s->f, 0);
//	 double y1 = gsl_vector_get (s->f, 1);
//	 double y2 = gsl_vector_get (s->f, 2);
//	 double y3 = gsl_vector_get (s->f, 3);
//	 std::cout << "the box width is :" << theBoxBounds.width() << std::endl;
	 
	 double t = gsl_vector_get(root,3);
	 double x0 = gsl_vector_get(root,0);
	 double x1 = gsl_vector_get(root,1);
	 double x2 = gsl_vector_get(root,2);	
	 
//	if (!in(t, Intervald(0.0,1.0))) status = GSL_CONTINUE;
	if(	status == GSL_SUCCESS) 
	{
		std::cout << "found critical point in the box:" << std::endl;
		theBoxBounds.print();

		std::cout << "t = " << t << ", at iter = " << iter << std::endl;
		TNTPoint critPoint(5);
		for(unsigned int j = 0; j <3; j++){
			critPoint[j] = gsl_vector_get (root, j);
		}
		
		critPoint[3] =  phi1 * t + (1 - t) * phi0; 
		critPoint[4] =	theta1 * t + (1 - t) * theta0; 
		std::cout << "phi = " << critPoint[3] << std::endl;
		std::cout << "theta = " << critPoint[4] << std::endl;
		
		gmVector3 cp(critPoint[0],critPoint[1],critPoint[2]);
		
		if(!SameRoot(critPoint, cpList)) 
			cpList.push_back(ClassifiedPoint(critPoint, singClassif->getSingularityType(cp)));
	}
#if 0	
	gsl_multiroot_fsolver_free (s);
#endif
	
#if 1	
	gsl_multiroot_fdfsolver_free (s);
#endif	
	gsl_vector_free (x);
	return status;
}

int
ContourCriticalPointsGSL::silh_kGkR_4DCP_df (const gsl_vector * xgsl, void *params, 
										gsl_matrix * J)
{
	//get parameters
	Implicit* F  = ((struct CPparams *) params)->F;
	double radius = ((struct CPparams *) params)->radius;
	double phi0 = ((struct CPparams *) params)->phi0; 
	double theta0 = ((struct CPparams *) params)->theta0;
	double phi1 = ((struct CPparams *) params)->phi1; 
	double theta1 = ((struct CPparams *) params)->theta1;
	
	//get gsl variables
	const double x0 = gsl_vector_get (xgsl, 0);
	const double x1 = gsl_vector_get (xgsl, 1);
	const double x2 = gsl_vector_get (xgsl, 2);
	const double x3 = gsl_vector_get (xgsl, 3);
	
	gmVector3 pos(x0,x1,x2);
	gmVector3 gradient(F->grad(pos));
	gmMatrix3 hessian(F->hess(pos));
	double theta = (1 - x3) * theta0 + x3 * theta1;
	double phi = (1 - x3) * phi0 + x3 * phi1;
	
	gmVector3 cameraPos( radius * sin(theta), -1 * radius * sin(phi) * cos(theta), -1 * radius * cos(phi) * cos(theta));	
	gmVector3 v(cameraPos - pos);
		
	double fx = F->Fx(pos);
	double fy = F->Fy(pos);
	double fz = F->Fz(pos);
	double fxx = F->Fxx(pos);
	double fxy = F->Fxy(pos);
	double fxz = F->Fxz(pos);
	double fyy = F->Fyy(pos);
	double fyz = F->Fyz(pos);
	double fzz = F->Fzz(pos);
	double fxxx = F->Fxxi(pos,0);
	double fxxy = F->Fxxi(pos,1);
	double fxyz = F->Fxyi(pos,2);
	double fyyy = F->Fyyi(pos,1);
	double fzzz = F->Fzzi(pos,2);
	double fxzz = F->Fxzi(pos,2);
	double fyyz = F->Fyyi(pos,2);
	double fxyy = F->Fxyi(pos,1);
	double fxxz = F->Fxxi(pos,2);
	double fyzz = F->Fyzi(pos,2);
	double denom = fx * fx + fy * fy + fz * fz;
	double x = x0;
	double y = x1;
	double z = x2;
	double cosPhi = cos(phi);
	double cosTheta = cos(theta);
	double sinPhi = sin(phi);
	double sinTheta = sin(theta);
	
//not an elegant avoidance of divide by 0 --ms 	
	if (denom == 0.0) 
	{
		double Jacobian[4][4] = { 
		{fx,fy,fz,0},
		{-fx + fxz*(-z - radius*cosPhi*cosTheta) + fxy*(-y - radius*cosTheta*sinPhi) + 
			fxx*(-x + radius*sinTheta),-fy + fyz*(-z - radius*cosPhi*cosTheta) + 
			fyy*(-y - radius*cosTheta*sinPhi) + fxy*(-x + radius*sinTheta),
			-fz + fzz*(-z - radius*cosPhi*cosTheta) + fyz*(-y - radius*cosTheta*sinPhi) + 
			fxz*(-x + radius*sinTheta),fx*radius*(-theta0 + theta1)*cosTheta + 
			fz*((-phi0 + phi1)*radius*cosTheta*sinPhi + radius*(-theta0 + theta1)*cosPhi*sinTheta) + 
			fy*(-((-phi0 + phi1)*radius*cosPhi*cosTheta) + radius*(-theta0 + theta1)*sinPhi*sinTheta)},
		{	HUGE, HUGE, HUGE },
		{-(fxz*(-z - radius*cosPhi*cosTheta)) - fxy*(-y - radius*cosTheta*sinPhi) - 
			fxx*(-x + radius*sinTheta) + (-x + radius*sinTheta)*
			(-fxx + fxxz*(-z - radius*cosPhi*cosTheta) + fxxy*(-y - radius*cosTheta*sinPhi) + 
			 fxxx*(-x + radius*sinTheta)) + (-y - radius*cosTheta*sinPhi)*
			(-fxy + fxyz*(-z - radius*cosPhi*cosTheta) + fxyy*(-y - radius*cosTheta*sinPhi) + 
			 fxxy*(-x + radius*sinTheta)) + (-z - radius*cosPhi*cosTheta)*
			(-fxz + fxzz*(-z - radius*cosPhi*cosTheta) + fxyz*(-y - radius*cosTheta*sinPhi) + 
			 fxxz*(-x + radius*sinTheta)),-(fyz*(-z - radius*cosPhi*cosTheta)) - 
			fyy*(-y - radius*cosTheta*sinPhi) - fxy*(-x + radius*sinTheta) + 
			(-x + radius*sinTheta)*(-fxy + fxyz*(-z - radius*cosPhi*cosTheta) + 
									fxyy*(-y - radius*cosTheta*sinPhi) + fxxy*(-x + radius*sinTheta)) + 
			(-y - radius*cosTheta*sinPhi)*(-fyy + fyyz*(-z - radius*cosPhi*cosTheta) + 
										   fyyy*(-y - radius*cosTheta*sinPhi) + fxyy*(-x + radius*sinTheta)) + 
			(-z - radius*cosPhi*cosTheta)*(-fyz + fyzz*(-z - radius*cosPhi*cosTheta) + 
										   fyyz*(-y - radius*cosTheta*sinPhi) + fxyz*(-x + radius*sinTheta)),
			-(fzz*(-z - radius*cosPhi*cosTheta)) - fyz*(-y - radius*cosTheta*sinPhi) - 
			fxz*(-x + radius*sinTheta) + (-x + radius*sinTheta)*
			(-fxz + fxzz*(-z - radius*cosPhi*cosTheta) + fxyz*(-y - radius*cosTheta*sinPhi) + 
			 fxxz*(-x + radius*sinTheta)) + (-y - radius*cosTheta*sinPhi)*
			(-fyz + fyzz*(-z - radius*cosPhi*cosTheta) + fyyz*(-y - radius*cosTheta*sinPhi) + 
			 fxyz*(-x + radius*sinTheta)) + (-z - radius*cosPhi*cosTheta)*
			(-fzz + fzzz*(-z - radius*cosPhi*cosTheta) + fyzz*(-y - radius*cosTheta*sinPhi) + 
			 fxzz*(-x + radius*sinTheta)),radius*(-theta0 + theta1)*cosTheta*
			(fxz*(-z - radius*cosPhi*cosTheta) + fxy*(-y - radius*cosTheta*sinPhi) + 
			 fxx*(-x + radius*sinTheta)) + (-((-phi0 + phi1)*radius*cosPhi*cosTheta) + 
											radius*(-theta0 + theta1)*sinPhi*sinTheta)*
			(fyz*(-z - radius*cosPhi*cosTheta) + fyy*(-y - radius*cosTheta*sinPhi) + 
			 fxy*(-x + radius*sinTheta)) + ((-phi0 + phi1)*radius*cosTheta*sinPhi + 
											radius*(-theta0 + theta1)*cosPhi*sinTheta)*
			(fzz*(-z - radius*cosPhi*cosTheta) + fyz*(-y - radius*cosTheta*sinPhi) + 
			 fxz*(-x + radius*sinTheta)) + (-x + radius*sinTheta)*
			(fxx*radius*(-theta0 + theta1)*cosTheta + 
			 fxz*((-phi0 + phi1)*radius*cosTheta*sinPhi + radius*(-theta0 + theta1)*cosPhi*sinTheta) + 
			 fxy*(-((-phi0 + phi1)*radius*cosPhi*cosTheta) + radius*(-theta0 + theta1)*sinPhi*sinTheta)) + 
			(-y - radius*cosTheta*sinPhi)*(fxy*radius*(-theta0 + theta1)*cosTheta + 
										   fyz*((-phi0 + phi1)*radius*cosTheta*sinPhi + radius*(-theta0 + theta1)*cosPhi*sinTheta) + 
										   fyy*(-((-phi0 + phi1)*radius*cosPhi*cosTheta) + radius*(-theta0 + theta1)*sinPhi*sinTheta)) + 
			(-z - radius*cosPhi*cosTheta)*(fxz*radius*(-theta0 + theta1)*cosTheta + 
										   fzz*((-phi0 + phi1)*radius*cosTheta*sinPhi + radius*(-theta0 + theta1)*cosPhi*sinTheta) + 
										   fyz*(-((-phi0 + phi1)*radius*cosPhi*cosTheta) + radius*(-theta0 + theta1)*sinPhi*sinTheta))}};	
		
		for(unsigned int i = 0; i<4; i++)
			for(unsigned int j = 0; j<4; j++)
				gsl_matrix_set (J, i, j, Jacobian[i][j]);
		
		return GSL_SUCCESS;
		
	}

///////////////////////
///Mathematica output	
///////////////////////	
	
	
	double Jacobian[4][4] = {{fx,fy,fz,0},{-fx + fxz*(-z - radius*cosPhi*cosTheta) + fxy*(-y - radius*cosTheta*sinPhi) + 
								   fxx*(-x + radius*sinTheta),-fy + fyz*(-z - radius*cosPhi*cosTheta) + 
								   fyy*(-y - radius*cosTheta*sinPhi) + fxy*(-x + radius*sinTheta),
								   -fz + fzz*(-z - radius*cosPhi*cosTheta) + fyz*(-y - radius*cosTheta*sinPhi) + 
								   fxz*(-x + radius*sinTheta),fx*radius*(-theta0 + theta1)*cosTheta + 
								   fz*((-phi0 + phi1)*radius*cosTheta*sinPhi + radius*(-theta0 + theta1)*cosPhi*sinTheta) + 
		fy*(-((-phi0 + phi1)*radius*cosPhi*cosTheta) + radius*(-theta0 + theta1)*sinPhi*sinTheta)},
{(-2.*(2.*fx*fxx + 2.*fxy*fy + 2.*fxz*fz)*((-pow(fxy,2) + fxx*fyy)*pow(fz,2) + 
														pow(fy,2)*(-pow(fxz,2) + fxx*fzz) + pow(fx,2)*(-pow(fyz,2) + fyy*fzz) - 
														2.*(fy*(-(fxy*fxz) + fxx*fyz)*fz + fx*(fxz*fyy - fxy*fyz)*fz + fx*fy*(-(fxz*fyz) + fxy*fzz))))/
				  pow(denom,3) + (2.*fxz*(-pow(fxy,2) + fxx*fyy)*fz + (-2.*fxxy*fxy + fxx*fxyy + fxxx*fyy)*pow(fz,2) + 
									2.*fxy*fy*(-pow(fxz,2) + fxx*fzz) + pow(fy,2)*(-2.*fxxz*fxz + fxx*fxzz + fxxx*fzz) + 
									pow(fx,2)*(fxzz*fyy - 2.*fxyz*fyz + fxyy*fzz) + 2.*fx*fxx*(-pow(fyz,2) + fyy*fzz) - 
									2.*(fxz*fy*(-(fxy*fxz) + fxx*fyz) + fx*fxz*(fxz*fyy - fxy*fyz) + fxy*(-(fxy*fxz) + fxx*fyz)*fz + 
										fy*(-(fxxz*fxy) + fxx*fxyz - fxxy*fxz + fxxx*fyz)*fz + 
										fx*(-(fxy*fxyz) + fxyy*fxz + fxxz*fyy - fxxy*fyz)*fz + fxx*(fxz*fyy - fxy*fyz)*fz + 
										fx*fy*(-(fxyz*fxz) + fxy*fxzz - fxxz*fyz + fxxy*fzz) + fx*fxy*(-(fxz*fyz) + fxy*fzz) + 
										fxx*fy*(-(fxz*fyz) + fxy*fzz)))/pow(denom,2),
				  (-2.*(2.*fx*fxy + 2.*fy*fyy + 2.*fyz*fz)*((-pow(fxy,2) + fxx*fyy)*pow(fz,2) + 
														pow(fy,2)*(-pow(fxz,2) + fxx*fzz) + pow(fx,2)*(-pow(fyz,2) + fyy*fzz) - 
														2.*(fy*(-(fxy*fxz) + fxx*fyz)*fz + fx*(fxz*fyy - fxy*fyz)*fz + fx*fy*(-(fxz*fyz) + fxy*fzz))))/
				  pow(denom,3) + (2.*(-pow(fxy,2) + fxx*fyy)*fyz*fz + (-2.*fxy*fxyy + fxxy*fyy + fxx*fyyy)*pow(fz,2) + 
									2.*fy*fyy*(-pow(fxz,2) + fxx*fzz) + pow(fy,2)*(-2.*fxyz*fxz + fxx*fyzz + fxxy*fzz) + 
									2.*fx*fxy*(-pow(fyz,2) + fyy*fzz) + pow(fx,2)*(-2.*fyyz*fyz + fyy*fyzz + fyyy*fzz) - 
									2.*(fy*fyz*(-(fxy*fxz) + fxx*fyz) + fx*fyz*(fxz*fyy - fxy*fyz) + fyy*(-(fxy*fxz) + fxx*fyz)*fz + 
										fy*(-(fxy*fxyz) - fxyy*fxz + fxx*fyyz + fxxy*fyz)*fz + fxy*(fxz*fyy - fxy*fyz)*fz + 
										fx*(fxyz*fyy + fxz*fyyy - fxy*fyyz - fxyy*fyz)*fz + fxy*fy*(-(fxz*fyz) + fxy*fzz) + 
										fx*fyy*(-(fxz*fyz) + fxy*fzz) + fx*fy*(-(fxz*fyyz) - fxyz*fyz + fxy*fyzz + fxyy*fzz)))/pow(denom,2),
				  (-2.*(2.*fx*fxz + 2.*fy*fyz + 2.*fz*fzz)*((-pow(fxy,2) + fxx*fyy)*pow(fz,2) + 
														pow(fy,2)*(-pow(fxz,2) + fxx*fzz) + pow(fx,2)*(-pow(fyz,2) + fyy*fzz) - 
														2.*(fy*(-(fxy*fxz) + fxx*fyz)*fz + fx*(fxz*fyy - fxy*fyz)*fz + fx*fy*(-(fxz*fyz) + fxy*fzz))))/
				  pow(denom,3) + ((-2.*fxy*fxyz + fxxz*fyy + fxx*fyyz)*pow(fz,2) + 2.*(-pow(fxy,2) + fxx*fyy)*fz*fzz + 
									2.*fy*fyz*(-pow(fxz,2) + fxx*fzz) + 2.*fx*fxz*(-pow(fyz,2) + fyy*fzz) + 
									pow(fy,2)*(-2.*fxz*fxzz + fxxz*fzz + fxx*fzzz) + pow(fx,2)*(-2.*fyz*fyzz + fyyz*fzz + fyy*fzzz) - 
									2.*(fyz*(-(fxy*fxz) + fxx*fyz)*fz + fxz*(fxz*fyy - fxy*fyz)*fz + 
										fy*(-(fxyz*fxz) - fxy*fxzz + fxxz*fyz + fxx*fyzz)*fz + 
										fx*(fxzz*fyy + fxz*fyyz - fxyz*fyz - fxy*fyzz)*fz + fy*(-(fxy*fxz) + fxx*fyz)*fzz + 
										fx*(fxz*fyy - fxy*fyz)*fzz + fxz*fy*(-(fxz*fyz) + fxy*fzz) + fx*fyz*(-(fxz*fyz) + fxy*fzz) + 
										fx*fy*(-(fxzz*fyz) - fxz*fyzz + fxyz*fzz + fxy*fzzz)))/pow(denom,2),0},
{-(fxz*(-z - radius*cosPhi*cosTheta)) - fxy*(-y - radius*cosTheta*sinPhi) - 
				  fxx*(-x + radius*sinTheta) + (-x + radius*sinTheta)*
				  (-fxx + fxxz*(-z - radius*cosPhi*cosTheta) + fxxy*(-y - radius*cosTheta*sinPhi) + 
				   fxxx*(-x + radius*sinTheta)) + (-y - radius*cosTheta*sinPhi)*
				  (-fxy + fxyz*(-z - radius*cosPhi*cosTheta) + fxyy*(-y - radius*cosTheta*sinPhi) + 
				   fxxy*(-x + radius*sinTheta)) + (-z - radius*cosPhi*cosTheta)*
				  (-fxz + fxzz*(-z - radius*cosPhi*cosTheta) + fxyz*(-y - radius*cosTheta*sinPhi) + 
				   fxxz*(-x + radius*sinTheta)),-(fyz*(-z - radius*cosPhi*cosTheta)) - 
				  fyy*(-y - radius*cosTheta*sinPhi) - fxy*(-x + radius*sinTheta) + 
				  (-x + radius*sinTheta)*(-fxy + fxyz*(-z - radius*cosPhi*cosTheta) + 
											fxyy*(-y - radius*cosTheta*sinPhi) + fxxy*(-x + radius*sinTheta)) + 
				  (-y - radius*cosTheta*sinPhi)*(-fyy + fyyz*(-z - radius*cosPhi*cosTheta) + 
													 fyyy*(-y - radius*cosTheta*sinPhi) + fxyy*(-x + radius*sinTheta)) + 
				  (-z - radius*cosPhi*cosTheta)*(-fyz + fyzz*(-z - radius*cosPhi*cosTheta) + 
													 fyyz*(-y - radius*cosTheta*sinPhi) + fxyz*(-x + radius*sinTheta)),
				  -(fzz*(-z - radius*cosPhi*cosTheta)) - fyz*(-y - radius*cosTheta*sinPhi) - 
				  fxz*(-x + radius*sinTheta) + (-x + radius*sinTheta)*
				  (-fxz + fxzz*(-z - radius*cosPhi*cosTheta) + fxyz*(-y - radius*cosTheta*sinPhi) + 
				   fxxz*(-x + radius*sinTheta)) + (-y - radius*cosTheta*sinPhi)*
				  (-fyz + fyzz*(-z - radius*cosPhi*cosTheta) + fyyz*(-y - radius*cosTheta*sinPhi) + 
				   fxyz*(-x + radius*sinTheta)) + (-z - radius*cosPhi*cosTheta)*
				  (-fzz + fzzz*(-z - radius*cosPhi*cosTheta) + fyzz*(-y - radius*cosTheta*sinPhi) + 
				   fxzz*(-x + radius*sinTheta)),radius*(-theta0 + theta1)*cosTheta*
				  (fxz*(-z - radius*cosPhi*cosTheta) + fxy*(-y - radius*cosTheta*sinPhi) + 
				   fxx*(-x + radius*sinTheta)) + (-((-phi0 + phi1)*radius*cosPhi*cosTheta) + 
													radius*(-theta0 + theta1)*sinPhi*sinTheta)*
				  (fyz*(-z - radius*cosPhi*cosTheta) + fyy*(-y - radius*cosTheta*sinPhi) + 
				   fxy*(-x + radius*sinTheta)) + ((-phi0 + phi1)*radius*cosTheta*sinPhi + 
													radius*(-theta0 + theta1)*cosPhi*sinTheta)*
				  (fzz*(-z - radius*cosPhi*cosTheta) + fyz*(-y - radius*cosTheta*sinPhi) + 
				   fxz*(-x + radius*sinTheta)) + (-x + radius*sinTheta)*
				  (fxx*radius*(-theta0 + theta1)*cosTheta + 
				   fxz*((-phi0 + phi1)*radius*cosTheta*sinPhi + radius*(-theta0 + theta1)*cosPhi*sinTheta) + 
				   fxy*(-((-phi0 + phi1)*radius*cosPhi*cosTheta) + radius*(-theta0 + theta1)*sinPhi*sinTheta)) + 
				  (-y - radius*cosTheta*sinPhi)*(fxy*radius*(-theta0 + theta1)*cosTheta + 
													 fyz*((-phi0 + phi1)*radius*cosTheta*sinPhi + radius*(-theta0 + theta1)*cosPhi*sinTheta) + 
													 fyy*(-((-phi0 + phi1)*radius*cosPhi*cosTheta) + radius*(-theta0 + theta1)*sinPhi*sinTheta)) + 
				  (-z - radius*cosPhi*cosTheta)*(fxz*radius*(-theta0 + theta1)*cosTheta + 
													 fzz*((-phi0 + phi1)*radius*cosTheta*sinPhi + radius*(-theta0 + theta1)*cosPhi*sinTheta) + 
													 fyz*(-((-phi0 + phi1)*radius*cosPhi*cosTheta) + radius*(-theta0 + theta1)*sinPhi*sinTheta))}};	
	
//////////////////////	
	
		
	for(unsigned int i = 0; i<4; i++)
		for(unsigned int j = 0; j<4; j++)
			gsl_matrix_set (J, i, j, Jacobian[i][j]);

	return GSL_SUCCESS;
}

int ContourCriticalPointsGSL::silh_kGkR_4DCP_fdf (const gsl_vector * x, void *params, gsl_vector * f, gsl_matrix * J)
{
	silh_kGkR_4DCP_f (x, params, f);
	silh_kGkR_4DCP_df (x, params, J);
	
	return GSL_SUCCESS;
}

double ContourCriticalPointsGSL::customWidth(const Box<double> & X)
{
	double maxwidth = -1.0;
	int dim, idx;
	
	for (idx = 0; idx < 3; idx++) 
		if (X[idx].width() > maxwidth) 
		{
			dim = idx;
			maxwidth = X[idx].width();
		}
			
	return 	maxwidth;	
}


// Function and Jacobian matrix for the Krawczyk contractor
INTERVAL_VECTOR fK(const INTERVAL_VECTOR& X, void* params)
{
	INTERVAL_VECTOR F(4);

	Implicit* Func  = ((struct CPparams *) params)->F;
	double radius = ((struct CPparams *) params)->radius;
	double phi0 = ((struct CPparams *) params)->phi0; 
	double theta0 = ((struct CPparams *) params)->theta0;
	double phi1 = ((struct CPparams *) params)->phi1; 
	double theta1 = ((struct CPparams *) params)->theta1;	
	
	Box3d Pos(X[0],X[1],X[2]);
	Box3d Gradient(Func->grad(Pos));
	Intervald theta = (1 - X[3]) * theta0 + X[3] * theta1;
	Intervald phi = (1 - X[3]) * phi0 + X[3] * phi1;

	Box3d CameraPos( radius * sin(theta), (-1) * radius * sin(phi) * cos(theta), -1 * radius * cos(phi) * cos(theta));
	
	Box3d V(CameraPos - Pos);
	Box3d P = Func->hess(Pos)*V;
	Intervald Kr = dot(P,V);
		
	F[0] = Func->proc(Pos);
	F[1] = dot(Gradient,V);
	F[2] = Func->gaussianCurvature(Pos); 
	F[3] = Kr;
	return F;
}

INTERVAL_VECTOR f5K(const INTERVAL_VECTOR& X, void* params)
{
	INTERVAL_VECTOR F(5);	
	Implicit* Func  = ((struct CPparams *) params)->F;
	double radius = ((struct CPparams *) params)->radius;
	
	Box3d Pos(X[0],X[1],X[2]);
	Box3d Gradient(Func->grad(Pos));
	Box3d CameraPos( radius * sin(X[4]), (-1) * radius * sin(X[3]) * cos(X[4]), (-1) * radius * cos(X[3]) * cos(X[4]));
	
	Box3d V(CameraPos - Pos);
	Box3d W(Func->hess(Pos) * V);
	Box3d G(Gradient[1]*W[2]- Gradient[2]*W[1],Gradient[2]*W[0]-
			Gradient[0]*W[2],Gradient[0]*W[1]- Gradient[1]*W[0]);
	
	F[0] = Func->proc(Pos);
	F[1] = dot(Gradient,V);
	F[2] = G[0];
	F[3] = G[1];
	F[4] = G[2];
	
	return F;
	
}

INTERVAL_MATRIX Jf5K(const INTERVAL_VECTOR& X, void* params)
{
	Implicit* F  = ((struct CPparams *) params)->F;
	double radius = ((struct CPparams *) params)->radius;
	
	Box3d Pos(X[0],X[1],X[2]);
	Box3d Gradient(F->grad(Pos));
	Box3d CameraPos(radius * sin(X[4]), (-1) * radius * sin(X[3]) * cos(X[4]),
			radius * cos(X[3]) * cos(X[4]));
	Box3d V(CameraPos - Pos);
	Box3d W(F->hess(Pos) * V);	
	IMatrix Hessian(F->hess(Pos));
	
	Box3d DX3CameraPos(0.0, (-1) * radius * cos(X[3]) * cos(X[4]),
			    radius * sin(X[3]) * cos(X[4]));
	Box3d DX4CameraPos(radius * cos(X[4]), radius * sin(X[3]) * sin(X[4]),
			    radius * cos(X[3]) * sin(X[4]));
	
	INTERVAL_MATRIX Jacobian(5,5); 
	
	for(unsigned int j = 0; j<3; j++) 
		Jacobian[0][j] = Gradient[j]; 
	
	Jacobian[0][3] = Jacobian[0][4] = 0.0;
	
	for(unsigned int j = 0; j<3; j++) 
		Jacobian[1][j] = W[j] - Gradient[j]; 
	
	Jacobian[1][3] = dot(Gradient,DX3CameraPos);	
	Jacobian[1][4] = dot(Gradient,DX4CameraPos);
	
	for(unsigned int i = 2; i<5; i++)
		for(unsigned int m = 0; m<3; m++)
			for(unsigned int j = 0; j<3; j++)
				for(unsigned int k = 0; k<3; k++)
					for(unsigned int l = 0; l<3; l++)
						Jacobian(i+1,m+1) += 
								epsilon[i-2][j][k] * Hessian(m+1,j+1) * Hessian(k+1,l+1) * V[l] + 
								epsilon[i-2][j][k] * Gradient[j] * F->D3F(Pos,m,k,l) * V[l] - epsilon[i-2][j][k] * Gradient[j] * Hessian(k+1,m+1); 
	
	for(unsigned int i = 2; i<5; i++)
		for(unsigned int j = 0; j<3; j++)
			for(unsigned int k = 0; k<3; k++)
				for(unsigned int l = 0; l<3; l++){
		Jacobian(i+1,4) += epsilon[i-2][j][k] * Gradient[j] * Hessian(k+1,l+1) * DX3CameraPos[l];
		Jacobian(i+1,5) += epsilon[i-2][j][k] * Gradient[j] * Hessian(k+1,l+1) * DX4CameraPos[l];
				}
					
	return Jacobian;
}


INTERVAL_MATRIX JfK(const INTERVAL_VECTOR& X, void* params)
{
	
//get parameters
	Implicit* F  = ((struct CPparams *) params)->F;
	double radius = ((struct CPparams *) params)->radius;
	double phi0 = ((struct CPparams *) params)->phi0; 
	double theta0 = ((struct CPparams *) params)->theta0;
	double phi1 = ((struct CPparams *) params)->phi1; 
	double theta1 = ((struct CPparams *) params)->theta1;	
	
	gmVector3 CameraPos0( radius * sin(theta0), (-1) * radius * sin(phi0) * cos(theta0), -1 * radius * cos(phi0) * cos(theta0));
	gmVector3 CameraPos1( radius * sin(theta1), (-1) * radius * sin(phi1) * cos(theta1), -1 * radius * cos(phi1) * cos(theta1));
	
	Box3d pos(X[0],X[1],X[2]);
	Intervald cp0 = (1 - X[3]) * CameraPos0[0] + X[3] * CameraPos1[0];
	Intervald cp1 = (1 - X[3]) * CameraPos0[1] + X[3] * CameraPos1[1];
	Intervald cp2 = (1 - X[3]) * CameraPos0[2] + X[3] * CameraPos1[2];
	
	Box3d cameraPos( cp0,cp1,cp2);
	Box3d gradient(F->grad(pos));
	IMatrix hessian(F->hess(pos));
	Box3d v(cameraPos - pos);
	
	Intervald fx = F->Fx(pos);
	Intervald fy = F->Fy(pos);
	Intervald fz = F->Fz(pos);
	Intervald fxx = F->Fxx(pos);
	Intervald fxy = F->Fxy(pos);
	Intervald fxz = F->Fxz(pos);
	Intervald fyy = F->Fyy(pos);
	Intervald fyz = F->Fyz(pos);
	Intervald fzz = F->Fzz(pos);
	Intervald fxxx = F->Fxxi(pos,0);
	Intervald fxxy = F->Fxxi(pos,1);
	Intervald fxyz = F->Fxyi(pos,2);
	Intervald fyyy = F->Fyyi(pos,1);
	Intervald fzzz = F->Fzzi(pos,2);
	Intervald fxzz = F->Fxzi(pos,2);
	Intervald fyyz = F->Fyyi(pos,2);
	Intervald fxyy = F->Fxyi(pos,1);
	Intervald fxxz = F->Fxxi(pos,2);
	Intervald fyzz = F->Fyzi(pos,2);
	Intervald denom = fx.squared() + fy.squared() + fz.squared();
	Intervald x = X[0];
	Intervald y = X[1];
	Intervald z = X[2];
	Intervald vx = v[0];
	Intervald vy = v[1];
	Intervald vz = v[2];
	double cx0 = CameraPos0[0];
	double cy0 = CameraPos0[1];
	double cz0 = CameraPos0[2];
	double cx1 = CameraPos1[0];
	double cy1 = CameraPos1[1];
	double cz1 = CameraPos1[2];

	INTERVAL_MATRIX Jacobian(4,4); 

//Mathematica output

	Jacobian[0][0] = fx;
	Jacobian[0][1] = fy;
	Jacobian[0][2] = fz;
	Jacobian[0][3] = 0;

	Jacobian[1][0] = -fx + fxx*vx + fxy*vy + fxz*vz;
	Jacobian[1][1] = -fy + fxy*vx + fyy*vy + fyz*vz;
	Jacobian[1][2] = -fz + fxz*vx + fyz*vy + fzz*vz;
	Jacobian[1][3] = (-cx0 + cx1)*fx + (-cy0 + cy1)*fy + 
				(-cz0 + cz1)*fz;
	Jacobian[2][0] = (-2*(2*fx*fxx + 2*fxy*fy + 2*fxz*fz)*
				   ((-fxy.squared() + fxx*fyy)*fz.squared() + 
					fy.squared()*(-fxz.squared() + fxx*fzz) + 
					fx.squared()*(-fyz.squared() + fyy*fzz) - 
					2.*(fy*(-(fxy*fxz) + fxx*fyz)*fz + 
						fx*(fxz*fyy - fxy*fyz)*fz + 
						fx*fy*(-(fxz*fyz) + fxy*fzz))))/
				  denom.pow(3) + 
				  (2*fxz*(-fxy.squared() + fxx*fyy)*fz + 
				   (-2*fxxy*fxy + fxx*fxyy + fxxx*fyy)*
				   fz.squared() + 
				   2*fxy*fy*(-fxz.squared() + fxx*fzz) + 
				   fy.squared()*(-2*fxxz*fxz + fxx*fxzz + 
								fxxx*fzz) + 
				   fx.squared()*(fxzz*fyy - 2*fxyz*fyz + 
								fxyy*fzz) + 
				   2*fx*fxx*(-fyz.squared() + fyy*fzz) - 
				   2.*(fxz*fy*(-(fxy*fxz) + fxx*fyz) + 
					   fx*fxz*(fxz*fyy - fxy*fyz) + 
					   fxy*(-(fxy*fxz) + fxx*fyz)*fz + 
					   fy*(-(fxxz*fxy) + fxx*fxyz - fxxy*fxz + 
						   fxxx*fyz)*fz + 
					   fx*(-(fxy*fxyz) + fxyy*fxz + fxxz*fyy - 
						   fxxy*fyz)*fz + 
					   fxx*(fxz*fyy - fxy*fyz)*fz + 
					   fx*fy*(-(fxyz*fxz) + fxy*fxzz - 
							  fxxz*fyz + fxxy*fzz) + 
					   fx*fxy*(-(fxz*fyz) + fxy*fzz) + 
					   fxx*fy*(-(fxz*fyz) + fxy*fzz)))/
		denom.squared();
		
		Jacobian[2][1] = (-2*(2*fx*fxy + 2*fy*fyy + 2*fyz*fz)*
				   ((-fxy.squared() + fxx*fyy)*fz.squared() + 
					fy.squared()*(-fxz.squared() + fxx*fzz) + 
					fx.squared()*(-fyz.squared() + fyy*fzz) - 
					2.*(fy*(-(fxy*fxz) + fxx*fyz)*fz + 
						fx*(fxz*fyy - fxy*fyz)*fz + 
						fx*fy*(-(fxz*fyz) + fxy*fzz))))/
				  denom.pow(3) + 
				  (2*(-fxy.squared() + fxx*fyy)*fyz*fz + 
				   (-2*fxy*fxyy + fxxy*fyy + fxx*fyyy)*
				   fz.squared() + 
				   2*fy*fyy*(-fxz.squared() + fxx*fzz) + 
				   fy.squared()*(-2*fxyz*fxz + fxx*fyzz + 
								fxxy*fzz) + 
				   2*fx*fxy*(-fyz.squared() + fyy*fzz) + 
				   fx.squared()*(-2*fyyz*fyz + fyy*fyzz + 
								fyyy*fzz) - 
				   2.*(fy*fyz*(-(fxy*fxz) + fxx*fyz) + 
					   fx*fyz*(fxz*fyy - fxy*fyz) + 
					   fyy*(-(fxy*fxz) + fxx*fyz)*fz + 
					   fy*(-(fxy*fxyz) - fxyy*fxz + fxx*fyyz + 
						   fxxy*fyz)*fz + 
					   fxy*(fxz*fyy - fxy*fyz)*fz + 
					   fx*(fxyz*fyy + fxz*fyyy - fxy*fyyz - 
						   fxyy*fyz)*fz + 
					   fxy*fy*(-(fxz*fyz) + fxy*fzz) + 
					   fx*fyy*(-(fxz*fyz) + fxy*fzz) + 
					   fx*fy*(-(fxz*fyyz) - fxyz*fyz + 
							  fxy*fyzz + fxyy*fzz)))/denom.squared();
		
	Jacobian[2][2] = (-2*(2*fx*fxz + 2*fy*fyz + 2*fz*fzz)*
					((-fxy.squared() + fxx*fyy)*fz.squared() + 
					 fy.squared()*(-fxz.squared() + fxx*fzz) + 
					 fx.squared()*(-fyz.squared() + fyy*fzz) - 
					 2.*(fy*(-(fxy*fxz) + fxx*fyz)*fz + 
						 fx*(fxz*fyy - fxy*fyz)*fz + 
						 fx*fy*(-(fxz*fyz) + fxy*fzz))))/
				  denom.pow(3) + 
				  ((-2*fxy*fxyz + fxxz*fyy + fxx*fyyz)*
				   fz.squared() + 
				   2*(-fxy.squared() + fxx*fyy)*fz*fzz + 
				   2*fy*fyz*(-fxz.squared() + fxx*fzz) + 
				   2*fx*fxz*(-fyz.squared() + fyy*fzz) + 
				   fy.squared()*(-2*fxz*fxzz + fxxz*fzz + 
								fxx*fzzz) + 
				   fx.squared()*(-2*fyz*fyzz + fyyz*fzz + 
								fyy*fzzz) - 
				   2.*(fyz*(-(fxy*fxz) + fxx*fyz)*fz + 
					   fxz*(fxz*fyy - fxy*fyz)*fz + 
					   fy*(-(fxyz*fxz) - fxy*fxzz + fxxz*fyz + 
						   fxx*fyzz)*fz + 
					   fx*(fxzz*fyy + fxz*fyyz - fxyz*fyz - 
						   fxy*fyzz)*fz + 
					   fy*(-(fxy*fxz) + fxx*fyz)*fzz + 
					   fx*(fxz*fyy - fxy*fyz)*fzz + 
					   fxz*fy*(-(fxz*fyz) + fxy*fzz) + 
					   fx*fyz*(-(fxz*fyz) + fxy*fzz) + 
					   fx*fy*(-(fxzz*fyz) - fxz*fyzz + 
							  fxyz*fzz + fxy*fzzz)))/denom.squared();
		Jacobian[2][3] = 0;

		Jacobian[3][0] = -(fxx*vx) - fxy*vy - fxz*vz + 
						   vx*(-fxx + fxxx*vx + fxxy*vy + fxxz*vz) + 
						   vy*(-fxy + fxxy*vx + fxyy*vy + fxyz*vz) + 
			vz*(-fxz + fxxz*vx + fxyz*vy + fxzz*vz);
			
		Jacobian[3][1] = -(fxy*vx) - fyy*vy - fyz*vz + 
						   vx*(-fxy + fxxy*vx + fxyy*vy + fxyz*vz) + 
						   vy*(-fyy + fxyy*vx + fyyy*vy + fyyz*vz) + 
			vz*(-fyz + fxyz*vx + fyyz*vy + fyzz*vz);
		
		Jacobian[3][2] = -(fxz*vx) - fyz*vy - fzz*vz + 
						   vx*(-fxz + fxxz*vx + fxyz*vy + fxzz*vz) + 
						   vy*(-fyz + fxyz*vx + fyyz*vy + fyzz*vz) + 
			vz*(-fzz + fxzz*vx + fyzz*vy + fzzz*vz);
		
		Jacobian[3][3] = ((-cx0 + cx1)*fxx + (-cy0 + cy1)*fxy + 
							(-cz0 + cz1)*fxz)*vx + 
						   ((-cx0 + cx1)*fxy + (-cy0 + cy1)*fyy + 
							(-cz0 + cz1)*fyz)*vy + 
						   ((-cx0 + cx1)*fxz + (-cy0 + cy1)*fyz + 
							(-cz0 + cz1)*fzz)*vz + 
						   (-cx0 + cx1)*(fxx*vx + fxy*vy + fxz*vz) + 
						   (-cy0 + cy1)*(fxy*vx + fyy*vy + fyz*vz) + 
							(-cz0 + cz1)*(fxz*vx + fyz*vy + fzz*vz);	
		return Jacobian;
}


