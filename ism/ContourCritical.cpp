/**
 * Implementation of a class used to find the critical points of a feature contour
 * using pure interval arithmetic methods.
 * @file ContourCritical.cpp
 * @date 01/13/2006
 * @author Matei N. Stroila
 * @remarks
 */

//by default, use a linear interpolation of the camera position
//if want to use an exact function in the solver, uncomment this:
//#define EXACT_CPS

#include "ContourCritical.h"
#include "Surface/gmTNTconvert.h"
#include "GaussJordan.h"
#include "Surface/Implicit/Operator/RotatedImplicit.h"
#include "SearchCritical.h"
#include "SingularityClassification.h"

/** 
* Default constructor.
*/
ContourCritical::ContourCritical()
{
}

ContourCritical::ContourCritical(CPparams *_params, double _NewtonTolerance, bool _useNewton)
{
	params = _params;
	NewtonTolerance = _NewtonTolerance;
	fUseNewton = _useNewton;
	WidthTolerance = 0.001;
	SameTolerance = WidthTolerance * 10.0;
	pRF = new RotatedImplicit();
	pRF->setF(params->F);
	params->pRF = pRF;
	params->pRinit = &Rinit;
	params->pRfinal = &Rfinal;
	singClassif = new SingularityClassification(pRF);
	
}

/**
* Overridden NewtonEquation to solve A(Y-Xc) = b for contour critical points.
 */
bool ContourCritical::NewtonEquation(Box<double> &X, IMatrix &A, Box<double> &b)
{
	//get the value of the function at the center of the interval
	TNTPoint fc(4);
	TNTPoint x = X.center();

#ifdef	EXACT_CPS	
	silh_kGkR_4DCP_f (x, fc);
#endif
#ifndef	EXACT_CPS		
	silh_interp_f (x, fc);
#endif
	
	//set the Jacobian
	IMatrix Jacobian(4,4); 
#ifdef	EXACT_CPS		
	silh_kGkR_4DCP_DF(X, Jacobian);
#endif
#ifndef	EXACT_CPS	
	silh_interp_DF(X, Jacobian);
#endif	
	
	TNT::Matrix<double> Vc = Jacobian.center();
	TNT::Matrix<double> VcInv(Vc.num_rows(),Vc.num_cols());
	bool okay = GaussJordan::Invert(Vc,VcInv);
	
	if (okay)
    {
		A = IMatrix(VcInv) * Jacobian;
		b = (-VcInv) * fc;
    }
	else // Vc was singular - couldn't invert it - don't use it
    {
		A = Jacobian;
		b = -fc;
    }
	return true;
}

/**
* Overridden Newton solver method to stop checking the current box for
 * roots during the subdivision phase.  
 */
bool ContourCritical::NewtonSubdivisionBreak(Box<double> & X)
{
	bool retval = false;
	int i;
	Interval<double> Fi;
	for (i = 0 ; i < X.size(); i++){
			//choose an exact function or an approximation 
			#ifdef	EXACT_CPS
				silh_kGkR_4DCP_F(X,Fi,i);
			#endif
			#ifndef	EXACT_CPS
				silh_interp_F(X,Fi,i);
			#endif
				if (Fi.isNegative() || Fi.isPositive()){
					retval =  true;   // No root in the region - we can quit early
					break;
				}
	}
	return retval;
}

void ContourCritical::silh_kGkR_4DCP_f (const TNTPoint& x, TNTPoint &f)
{
	//get parameters
	Implicit* F  = params->F;
	double radius = params->radius;
	double phi0 = params->phi0; 
	double theta0 = params->theta0;
	double phi1 = params->phi1; 
	double theta1 = params->theta1;
		
	gmVector3 pos(x[0],x[1],x[2]);
	gmVector3 gradient(F->grad(pos));
	double theta = (1 - x[3]) * theta0 + x[3] * theta1;
	double phi = (1 - x[3]) * phi0 + x[3] * phi1;
	
	gmVector3 cameraPos( radius * sin(theta), -1 * radius * sin(phi) * cos(theta), -1 * radius * cos(phi) * cos(theta));
	
	gmVector3 v(cameraPos - pos);
	
	f[0] = F->proc(pos);
	f[1] = dot(gradient,v);
	f[2] = F->gaussianCurvature(pos); 
	f[3] = dot(F->hess(pos)*v,v);
	
}

void ContourCritical::silh_kGkR_4DCP_F(const Box<double>& X, Interval<double>& Fi, const int i)
{
	Implicit* Func  = params->F;
	double radius = params->radius;
	double phi0 = params->phi0; 
	double theta0 = params->theta0;
	double phi1 = params->phi1; 
	double theta1 = params->theta1;	
	
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


void ContourCritical::silh_kGkR_4DCP_DF(const Box<double> & X, IMatrix& Jacobian)
{
	
	//get parameters
	Implicit* F  = params->F;
	double radius = params->radius;
	double phi0 = params->phi0; 
	double theta0 = params->theta0;
	double phi1 = params->phi1; 
	double theta1 = params->theta1;
	
	Box3d pos(X[0],X[1],X[2]);
	Box3d gradient(F->grad(pos));
	IMatrix hessian(F->hess(pos));
	Intervald theta = (1 - X[3]) * theta0 + X[3] * theta1;
	Intervald phi = (1 - X[3]) * phi0 + X[3] * phi1;
	Box3d cameraPos( radius * sin(theta), -1 * radius * sin(phi) * cos(theta), -1 * radius * cos(phi) * cos(theta));	
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
	Intervald cosPhi = cos(phi);
	Intervald cosTheta = cos(theta);
	Intervald sinPhi = sin(phi);
	Intervald sinTheta = sin(theta);
	
	///////////////////////
	///Mathematica output	
	///////////////////////		
	
	Jacobian[0][0] = fx;
	Jacobian[0][1] = fy;
	Jacobian[0][2] = fz;
	Jacobian[0][3] = 0;
	
	Jacobian[1][0] = -fx + fxz*(-z - radius*cosPhi*cosTheta) + fxy*(-y - radius*cosTheta*sinPhi) + fxx*(-x + radius*sinTheta);
	Jacobian[1][1] = -fy + fyz*(-z - radius*cosPhi*cosTheta) + 
		fyy*(-y - radius*cosTheta*sinPhi) + fxy*(-x + radius*sinTheta);
	Jacobian[1][2] = -fz + fzz*(-z - radius*cosPhi*cosTheta) + fyz*(-y - radius*cosTheta*sinPhi) + fxz*(-x + radius*sinTheta);
	Jacobian[1][3] = fx*radius*(-theta0 + theta1)*cosTheta + fz*((-phi0 + phi1)*radius*cosTheta*sinPhi + radius*(-theta0 + theta1)*cosPhi*sinTheta) + fy*(-((-phi0 + phi1)*radius*cosPhi*cosTheta) + radius*(-theta0 + theta1)*sinPhi*sinTheta);
	
	Jacobian[2][0] = (-2.*(2.*fx*fxx + 2.*fxy*fy + 2.*fxz*fz)*((-fxy.squared() + fxx*fyy)*fz.squared() + fy.squared()*(-fxz.squared() + fxx*fzz) + fx.squared()*(-fyz.squared() + fyy*fzz) - 2.*(fy*(-(fxy*fxz) + fxx*fyz)*fz + fx*(fxz*fyy - fxy*fyz)*fz + fx*fy*(-(fxz*fyz) + fxy*fzz))))/
		denom.pow(3) + (2.*fxz*(-fxy.squared() + fxx*fyy)*fz + (-2.*fxxy*fxy + fxx*fxyy + fxxx*fyy)*fz.squared() + 2.*fxy*fy*(-fxz.squared() + fxx*fzz) + fy.squared()*(-2.*fxxz*fxz + fxx*fxzz + fxxx*fzz) + fx.squared()*(fxzz*fyy - 2.*fxyz*fyz + fxyy*fzz) + 2.*fx*fxx*(-fyz.squared() + fyy*fzz) - 2.*(fxz*fy*(-(fxy*fxz) + fxx*fyz) + fx*fxz*(fxz*fyy - fxy*fyz) + fxy*(-(fxy*fxz) + fxx*fyz)*fz + fy*(-(fxxz*fxy) + fxx*fxyz - fxxy*fxz + fxxx*fyz)*fz + 
  fx*(-(fxy*fxyz) + fxyy*fxz + fxxz*fyy - fxxy*fyz)*fz + fxx*(fxz*fyy - fxy*fyz)*fz + 
 fx*fy*(-(fxyz*fxz) + fxy*fxzz - fxxz*fyz + fxxy*fzz) + fx*fxy*(-(fxz*fyz) + fxy*fzz) + 
  fxx*fy*(-(fxz*fyz) + fxy*fzz)))/denom.squared();
	Jacobian[2][1] = (-2.*(2.*fx*fxy + 2.*fy*fyy + 2.*fyz*fz)*((-fxy.squared() + fxx*fyy)*fz.squared() +  fy.squared()*(-fxz.squared() + fxx*fzz) + fx.squared()*(-fyz.squared() + fyy*fzz) -  2.*(fy*(-(fxy*fxz) + fxx*fyz)*fz + fx*(fxz*fyy - fxy*fyz)*fz + fx*fy*(-(fxz*fyz) + fxy*fzz))))/
		denom.pow(3) + (2.*(-fxy.squared() + fxx*fyy)*fyz*fz + (-2.*fxy*fxyy + fxxy*fyy + fxx*fyyy)*fz.squared() + 2.*fy*fyy*(-fxz.squared() + fxx*fzz) + fy.squared()*(-2.*fxyz*fxz + fxx*fyzz + fxxy*fzz) +	2.*fx*fxy*(-fyz.squared() + fyy*fzz) + fx.squared()*(-2.*fyyz*fyz + fyy*fyzz + fyyy*fzz) -2.*(fy*fyz*(-(fxy*fxz) + fxx*fyz) + fx*fyz*(fxz*fyy - fxy*fyz) + fyy*(-(fxy*fxz) + fxx*fyz)*fz +	fy*(-(fxy*fxyz) - fxyy*fxz + fxx*fyyz + fxxy*fyz)*fz + fxy*(fxz*fyy - fxy*fyz)*fz + fx*(fxyz*fyy + fxz*fyyy - fxy*fyyz - fxyy*fyz)*fz + fxy*fy*(-(fxz*fyz) + fxy*fzz) + fx*fyy*(-(fxz*fyz) + fxy*fzz) + fx*fy*(-(fxz*fyyz) - fxyz*fyz + fxy*fyzz + fxyy*fzz)))/denom.squared();
	Jacobian[2][2] = (-2.*(2.*fx*fxz + 2.*fy*fyz + 2.*fz*fzz)*((-fxy.squared() + fxx*fyy)*fz.squared() + fy.squared()*(-fxz.squared() + fxx*fzz) + fx.squared()*(-fyz.squared() + fyy*fzz) - 2.*(fy*(-(fxy*fxz) + fxx*fyz)*fz + fx*(fxz*fyy - fxy*fyz)*fz + fx*fy*(-(fxz*fyz) + fxy*fzz))))/
		denom.pow(3) + ((-2.*fxy*fxyz + fxxz*fyy + fxx*fyyz)*fz.squared() + 2.*(-fxy.squared() + fxx*fyy)*fz*fzz + 2.*fy*fyz*(-fxz.squared() * fx.squared() + fxx*fzz) + 2.*fx*fxz*(-fyz.squared() + fyy*fzz) + 
						fy.squared()*(-2.*fxz*fxzz + fxxz*fzz + fxx*fzzz) + fx.squared()*(-2.*fyz*fyzz + fyyz*fzz + fyy*fzzz) - 2.*(fyz*(-(fxy*fxz) + fxx*fyz)*fz + fxz*(fxz*fyy - fxy*fyz)*fz + fy*(-(fxyz*fxz) - fxy*fxzz + fxxz*fyz + fxx*fyzz)*fz + fx*(fxzz*fyy + fxz*fyyz - fxyz*fyz - fxy*fyzz)*fz + fy*(-(fxy*fxz) + fxx*fyz)*fzz + fx*(fxz*fyy - fxy*fyz)*fzz + fxz*fy*(-(fxz*fyz) + fxy*fzz) + fx*fyz*(-(fxz*fyz) + fxy*fzz) + fx*fy*(-(fxzz*fyz) - fxz*fyzz + fxyz*fzz + fxy*fzzz)))/denom.squared();
	Jacobian[2][3] = 0;
	
	Jacobian[3][0] = -(fxz*(-z - radius*cosPhi*cosTheta)) - fxy*(-y - radius*cosTheta*sinPhi) - fxx*(-x + radius*sinTheta) + (-x + radius*sinTheta)*
		(-fxx + fxxz*(-z - radius*cosPhi*cosTheta) + fxxy*(-y - radius*cosTheta*sinPhi) + 
		 fxxx*(-x + radius*sinTheta)) + (-y - radius*cosTheta*sinPhi)*
		(-fxy + fxyz*(-z - radius*cosPhi*cosTheta) + fxyy*(-y - radius*cosTheta*sinPhi) + 
		 fxxy*(-x + radius*sinTheta)) + (-z - radius*cosPhi*cosTheta)*
		(-fxz + fxzz*(-z - radius*cosPhi*cosTheta) + fxyz*(-y - radius*cosTheta*sinPhi) + 
		 fxxz*(-x + radius*sinTheta));
	Jacobian[3][1] = -(fyz*(-z - radius*cosPhi*cosTheta)) - 
		fyy*(-y - radius*cosTheta*sinPhi) - fxy*(-x + radius*sinTheta) + 
		(-x + radius*sinTheta)*(-fxy + fxyz*(-z - radius*cosPhi*cosTheta) + fxyy*(-y - radius*cosTheta*sinPhi) + fxxy*(-x + radius*sinTheta)) + 
		(-y - radius*cosTheta*sinPhi)*(-fyy + fyyz*(-z - radius*cosPhi*cosTheta) + 
									   fyyy*(-y - radius*cosTheta*sinPhi) + fxyy*(-x + radius*sinTheta)) + 
		(-z - radius*cosPhi*cosTheta)*(-fyz + fyzz*(-z - radius*cosPhi*cosTheta) + 
									   fyyz*(-y - radius*cosTheta*sinPhi) + fxyz*(-x + radius*sinTheta));
	Jacobian[3][2] = -(fzz*(-z - radius*cosPhi*cosTheta)) - fyz*(-y - radius*cosTheta*sinPhi) - fxz*(-x + radius*sinTheta) + (-x + radius*sinTheta)*
		(-fxz + fxzz*(-z - radius*cosPhi*cosTheta) + fxyz*(-y - radius*cosTheta*sinPhi) + 
		 fxxz*(-x + radius*sinTheta)) + (-y - radius*cosTheta*sinPhi)*
		(-fyz + fyzz*(-z - radius*cosPhi*cosTheta) + fyyz*(-y - radius*cosTheta*sinPhi) + 
		 fxyz*(-x + radius*sinTheta)) + (-z - radius*cosPhi*cosTheta)*
		(-fzz + fzzz*(-z - radius*cosPhi*cosTheta) + fyzz*(-y - radius*cosTheta*sinPhi) + 
		 fxzz*(-x + radius*sinTheta));
	Jacobian[3][3] = radius*(-theta0 + theta1)*cosTheta*
		(fxz*(-z - radius*cosPhi*cosTheta) + fxy*(-y - radius*cosTheta*sinPhi) + 
		 fxx*(-x + radius*sinTheta)) + (-((-phi0 + phi1)*radius*cosPhi*cosTheta) +  radius*(-theta0 + theta1)*sinPhi*sinTheta)*
		(fyz*(-z - radius*cosPhi*cosTheta) + fyy*(-y - radius*cosTheta*sinPhi) + 
		 fxy*(-x + radius*sinTheta)) + ((-phi0 + phi1)*radius*cosTheta*sinPhi + radius*(-theta0 + theta1)*cosPhi*sinTheta)*
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
									   fyz*(-((-phi0 + phi1)*radius*cosPhi*cosTheta) + radius*(-theta0 + theta1)*sinPhi*sinTheta));	
	

}

void ContourCritical::silh_interp_f (const TNTPoint& x, TNTPoint &f)
{
	//get parameters
	Implicit* F  = params->F;
	double radius = params->radius;
	double phi0 = params->phi0; 
	double theta0 = params->theta0;
	double phi1 = params->phi1; 
	double theta1 = params->theta1;
	
	gmVector3 CameraPos0( radius * sin(theta0), (-1) * radius * sin(phi0) * cos(theta0), -1 * radius * cos(phi0) * cos(theta0));
	gmVector3 CameraPos1( radius * sin(theta1), (-1) * radius * sin(phi1) * cos(theta1), -1 * radius * cos(phi1) * cos(theta1));
	
	double cp0 = (1 - x[3]) * CameraPos0[0] + x[3] * CameraPos1[0];
	double cp1 = (1 - x[3]) * CameraPos0[1] + x[3] * CameraPos1[1];
	double cp2 = (1 - x[3]) * CameraPos0[2] + x[3] * CameraPos1[2];
	
	gmVector3 cameraPos( cp0,cp1,cp2);
	
	gmVector3 pos(x[0],x[1],x[2]);
	gmVector3 gradient(F->grad(pos));
	
	gmVector3 v(cameraPos - pos);
	
	f[0] = F->proc(pos);
	f[1] = dot(gradient,v);
	f[2] = F->gaussianCurvature(pos); 
	f[3] = dot(F->hess(pos)*v,v);
	
}


void ContourCritical::silh_interp_F(const Box<double>& X, Interval<double>& Fi, const int i)
{
	Implicit* Func  = params->F;
	double radius = params->radius;
	double phi0 = params->phi0; 
	double theta0 = params->theta0;
	double phi1 = params->phi1; 
	double theta1 = params->theta1;	
	
	gmVector3 CameraPos0( radius * sin(theta0), (-1) * radius * sin(phi0) * cos(theta0), -1 * radius * cos(phi0) * cos(theta0));
	gmVector3 CameraPos1( radius * sin(theta1), (-1) * radius * sin(phi1) * cos(theta1), -1 * radius * cos(phi1) * cos(theta1));
	
	Box3d Pos(X[0],X[1],X[2]);
	Box3d Gradient,V, P;

	Intervald cp0 = (1 - X[3]) * CameraPos0[0] + X[3] * CameraPos1[0];
	Intervald cp1 = (1 - X[3]) * CameraPos0[1] + X[3] * CameraPos1[1];
	Intervald cp2 = (1 - X[3]) * CameraPos0[2] + X[3] * CameraPos1[2];
	
	Box3d CameraPos( cp0,cp1,cp2);

	switch (i){
		case 0:
			Fi = Func->proc(Pos);
			break;
		case 1:
			Gradient= Func->grad(Pos);
			V = CameraPos - Pos;
			Fi = dot(Gradient,V);
		    break;
		case 2:
			Fi = Func->gaussianCurvature(Pos);
			break;
		case 3:
			V = CameraPos - Pos;
			P = Func->hess(Pos)*V;
			Fi = dot(P,V);
		    break;
		default:
			std::cout << "invalid i = " << i << std::endl;
	}
	
}

void ContourCritical::silh_interp_DF(const Box<double> & X, IMatrix& Jacobian)
{
	
	//get parameters
	Implicit* F  = params->F;
	double radius = params->radius;
	double phi0 = params->phi0; 
	double theta0 = params->theta0;
	double phi1 = params->phi1; 
	double theta1 = params->theta1;	
	
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
}	

int ContourCritical::search(Box<double> &X, ContourType type, ClassifiedPointList& cpList)
{
	PointList::const_iterator iter;
	fNewtonTolerance = NewtonTolerance;
	FindZeros(X);
	for(iter = fRoots.begin(); iter !=fRoots.end(); ++iter)
	{
		gmVector3 cp((*iter)[0],(*iter)[1],(*iter)[2]);
		cpList.push_back(ClassifiedPoint(*iter, singClassif->getSingularityType(cp)));
	}
	
//set 1 to check the newton solver on a simpler system (critical points of the implicit)
#if 0	 
	SearchCritical sc(params->F);
	sc.search();
	for(iter = sc.fRoots.begin(); iter !=sc.fRoots.end(); ++iter)
	{
		gmVector3 cp((*iter)[0],(*iter)[1],(*iter)[2]);
		cpList.push_back(ClassifiedPoint(*iter, singClassif->getSingularityType(cp)));
	}	
#endif

	
	return fRoots.size();
}

ContourCritical::~ContourCritical()
{
	delete pRF;
	delete singClassif;
}

