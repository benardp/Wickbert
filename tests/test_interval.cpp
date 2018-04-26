	//test the interval library
	
	#include <iostream>
	#include "Surface/Interval.h"
	#include "Surface/IMatrix.h"
	#include "ism/contractors.h"
	
	// Function and Jacobian matrix for the Krawczyk contractor
	INTERVAL_VECTOR fK(const INTERVAL_VECTOR& x, void* params)
	{
		INTERVAL_VECTOR r(2);
		r(1) = square(x(1)) - 1.0;
		r(2) = square(x(2)) - 1.0;
		
		return r;
	}
	
	INTERVAL_MATRIX JfK(const INTERVAL_VECTOR& x, void* params)
	{
		INTERVAL_MATRIX r(2,2);
		r(1,1) = 2 * x(1);    r(1,2) = 0;
		r(2,1) = 0;          r(2,2) = 2 * x(2);
		
		return r;
	}
	
	int FindZerosIA(Box<double> &X)
	{
		Box<double> X1,X2; /* The two halves if we divide X */
		double width = -1;      /* The width of the input box */
					
		int numberCP = 0;
	
//		X = KrawczykContractor(fK, JfK, X, NULL);
	
		int dim, idx;
	
		for (idx = 0; idx < 2; idx++) 
			if (X[idx].width() > width) 
		{
			dim = idx;
			width = X[idx].width();
		}
		
		if (width < 0.000001) {
			std::cout << "added critical point in the box:" << std::endl;
			X.print();
			return 1;
		}
		/* Unsure, so subdivide */	
		double maxwidth = -1.0;
		for (int idx = 0; idx < 2; idx++) 
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
		numberCP += FindZerosIA(X1);
		numberCP += FindZerosIA(X2);
		return numberCP;
	} 
	
	int main() {
	
		Intervald I(-1,1);
		std::cout << "(-1,1)^2= " << I.squared().low() << "," <<  I.squared().high() << std::endl;
		std::cout << "(-1,1) * (-1,1) = " << (I * I).low() << "," <<  (I * I).high() << std::endl;
		std::cout << "pow(I,2) = " << I.pow(2).low() << "," <<   I.pow(2).high() << std::endl;
		std::cout << "(-1,1) < 2.0 = " << (I < 2.0) << std::endl;
		std::cout << "(-1,1) ^ -1 " << I.pow(-1).low() << "," <<   I.pow(-1).high()  << std::endl;
	
		INTERVAL_MATRIX A(3,3);
		
		A(1,1) = INTERVAL(4,5);          A(1,2) = INTERVAL(-1,1);        A(1,3) = INTERVAL(1.5,2.5);
		A(2,1) = INTERVAL(-0.5,0.5); A(2,2) = INTERVAL(-7,-5);   A(2,3) = INTERVAL(1,2);
		A(3,1) = INTERVAL(-1.5,-0.5);   A(3,2) = INTERVAL(-0.7,-0.5); A(3,3) = INTERVAL(2,3);
		
		INTERVAL_VECTOR b(3);
		
		b(1) = INTERVAL(3,4);
		b(2) = INTERVAL(0,2);
		b(3) = INTERVAL(3,4);
		
		INTERVAL_VECTOR p(3);
		p(1) = INTERVAL(-10,10);
		p(2) = INTERVAL(-10,10);
		p(3) = INTERVAL(-10,10);
		
		std::cout << " A = " << A << std::endl;
		std::cout << " b = " << b << std::endl;
		std::cout << " p = " << p << std::endl;
		
		std::cout << "Solution of 'Ap-b=0' using Gauss elimination contractor\n" ;
		
		p = GaussElimination(A,b,p);
		
		std::cout.precision(5);
		std::cout << "p = " << p << std::endl;
		
		p(1) = INTERVAL(-10,10);
		p(2) = INTERVAL(-10,10);
		p(3) = INTERVAL(-10,10);
		
		std::cout << "Solution of 'Ap-b=0' using Gauss-Seidel contractor\n" ;
		
		std::cout.precision(5);
		for (int i=1;i<=100;i++)
		{
			p = GaussSeidelIteration(A,b,p);
			if ((i==1)||(i==2)||(i==5)||(i==10)||(i==20)||(i==100))
				std::cout << "Iteration " << i << " : p = " <<  p <<  std::endl;
		}
		
		std::cout << std::endl << "Press the 'return' key for contractors with preconditionning" << std::endl;
		std::cin.get();
		
		p(1) = INTERVAL(-10,10);
		p(2) = INTERVAL(-10,10);
		p(3) = INTERVAL(-10,10);
		
		std::cout << "Solution using Gauss elimination contractor with preconditionning\n" ;
		p = PrecGaussElimination(A,b,p);
		
		std::cout.precision(5);
		std::cout << "p = " << p << std::endl;
		
		
		p(1) = INTERVAL(-10,10);
		p(2) = INTERVAL(-10,10);
		p(3) = INTERVAL(-10,10);
		
		std::cout << "Solution using Gauss-Seidel contractor with preconditionning\n" ;
		
		std::cout.precision(5);
		for (int i=1;i<=100;i++)
		{
			p = PrecGaussSeidelIteration(A,b,p);
			if ((i==1)||(i==2)||(i==5)||(i==10)||(i==20)||(i==100))
			std::cout << "Iteration " << i << " : p = " <<  p << std::endl;
		}
		
		std::cout << std::endl << "Press the 'return' key for the Krawczyk contractor" << std::endl;
		std::cin.get();
		
		std::cout << "The problem is to solve the CSP H: (f(x)=0, with x in [x])" << std::endl;
		std::cout << "f1(x) = square(x(1)) - 4 * x(2)," << std::endl;
		std::cout << "f1(x) = square(x(2)) - 2 * x(1) + 4 * x(2)," << std::endl;
		std::cout << "[x] = [-0.1,0.1]x[-0.1,0.3]." << std::endl;
		
		INTERVAL_VECTOR x(2);
		x(1) = INTERVAL(-1.1,1.1);
		x(2) = INTERVAL(-1.1,1.3);
		
//		FindZerosIA(x);
		
		std::cin.get();

		std::cout << "box: " << x << std::endl;
		std::cout << "its center: " << x.center() << std::endl;	
		Box<double> xo(x.center());
		std::cout << "back to a box: " << xo << std::endl;	
		Box<double> diff(x - xo);
		std::cout << "x-xo: " << diff << std::endl;
		return 0;
	}
