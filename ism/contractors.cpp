/*
 *  contractors.cpp
 *  gps
 *
 *  Created by Matei Stroila on 2/27/06.
 *  
 *
 */


#include "contractors.h"
#include "GaussJordan.h"

// contractor based on Gauss elimination
INTERVAL_VECTOR GaussElimination(INTERVAL_MATRIX A,
        INTERVAL_VECTOR b, INTERVAL_VECTOR p)
{
	// Dimension of the problem
	int np = b.size();
	INTERVAL_VECTOR alpha(np);
	int i,j,k;
	for (i=1;i<=np-1;i++)
	{
		for (j=i+1;j<=np;j++)
		{
			alpha(j) = A(j,i)/A(i,i);
			b(j) = b(j)-alpha(j)*b(i);
		
			for (k=i+1;k<=np;k++)
				A(j,k) = A(j,k) - alpha(j)*A(i,k);
		}
	}
	
	for (i=np;i>=1;i--)
	{
		INTERVAL temp(0);
		for (j=i+1;j<=np;j++) temp += A(i,j)*p(j);
		
		temp = (b(i) - temp)/A(i,i);
		p(i) = p(i).intersection(temp);
	}
	
	return p;
}

// contractor based on Gauss-Seidel
INTERVAL_VECTOR GaussSeidelIteration(INTERVAL_MATRIX A,
        INTERVAL_VECTOR b, INTERVAL_VECTOR p)
{
	int np = b.size();
	INTERVAL_MATRIX AdiagInv(np,np),Aantidiag(A);
	INTERVAL_VECTOR prc(np,0);
	int i;
	
	for (i=1;i<=np;i++)
	{
		AdiagInv(i,i) = 1/A(i,i);
		Aantidiag(i,i) = 0;
	}
	
	INTERVAL_VECTOR temp = Aantidiag*p;
	for (i=0;i< np;i++)
	{
		temp[i] = b[i]-temp[i];
	}
	prc = AdiagInv * temp;
	
	for (i=0;i< np;i++)
	{
		prc[i] = prc[i].intersection(p[i]);
	}
	return (prc);
 }

// contractors using preconditioning
INTERVAL_VECTOR PrecGaussElimination(INTERVAL_MATRIX A,
        INTERVAL_VECTOR b, INTERVAL_VECTOR p)
{
	MATRIX AInv(A.num_rows(), A.num_cols());
	MATRIX Ac = A.center();
	GaussJordan::Invert(Ac, AInv);
	INTERVAL_MATRIX AInvI(AInv);
	INTERVAL_MATRIX Ap = AInvI*A;
	INTERVAL_VECTOR bp = AInvI*b;
	return GaussElimination(Ap,bp,p);
}


INTERVAL_VECTOR PrecGaussSeidelIteration(INTERVAL_MATRIX A,
        INTERVAL_VECTOR b, INTERVAL_VECTOR p)
{
   
	MATRIX AInv(A.num_rows(), A.num_cols());
	MATRIX Ac = A.center();
	if(!GaussJordan::Invert(Ac, AInv)) return GaussSeidelIteration(A,b,p);
	INTERVAL_MATRIX AInvI(AInv);
	INTERVAL_MATRIX Ap = AInvI*A;
	INTERVAL_VECTOR bp = AInvI*b;
  	return GaussSeidelIteration(Ap,bp,p);

}

// Krawczyk contractor
/*
INTERVAL_VECTOR KrawczykContractor(PIVF f, PIMF Jf, INTERVAL_VECTOR x, void* params)
{
	int np = x.size();
	INTERVAL_VECTOR x0(x.center());
	MATRIX M(np, np);
	MATRIX Jfc = Jf(x0, params).center();
	GaussJordan::Invert(Jfc, M);
	INTERVAL_MATRIX MI(M);
	INTERVAL_MATRIX Id(np,np);
	Id.MakeIdentity();
	INTERVAL_MATRIX Jpsi = Id - MI*Jf(x, params);
	TNT_INTERVAL_VECTOR r = x0-MI*f(x0, params)+Jpsi*(x-x0);
	INTERVAL_VECTOR retX(np);
	INTERVAL_VECTOR retX(np);
	for (int i=0;i< np;i++)
	{
		retX[i] = x[i].intersection(r[i]);
	}
	return retX;

}

*/
