/*
 *  contractors.h
 *  gps
 *
 *  Created by Matei Stroila on 2/27/06.
 *  
 *
 */
#ifndef contractors_h
#define contractors_h

#include "Surface/Box.h"
#include "Surface/IMatrix.h"

typedef Box<double> INTERVAL_VECTOR;
typedef IMatrix INTERVAL_MATRIX;
typedef Interval<double> INTERVAL;
typedef TNT::Matrix<double> MATRIX;
typedef TNT::Vector<INTERVAL> TNT_INTERVAL_VECTOR;

// defines "pointer to an INTERVAL_VECTOR function"
typedef INTERVAL_VECTOR (*PIVF)(const INTERVAL_VECTOR&, void* params);
// defines "pointer to an INTERVAL_MATRIX function"
typedef INTERVAL_MATRIX (*PIMF)(const INTERVAL_VECTOR&, void* params);

// contractor based on Gauss elimination
INTERVAL_VECTOR GaussElimination(INTERVAL_MATRIX A,
        INTERVAL_VECTOR b, INTERVAL_VECTOR p);

// contractor based on Gauss-Seidel elimination
INTERVAL_VECTOR GaussSeidelIteration(INTERVAL_MATRIX A,
        INTERVAL_VECTOR b, INTERVAL_VECTOR p);

// contractors using preconditioning
INTERVAL_VECTOR PrecGaussElimination(INTERVAL_MATRIX A,
        INTERVAL_VECTOR b, INTERVAL_VECTOR p);

INTERVAL_VECTOR PrecGaussSeidelIteration(INTERVAL_MATRIX A,
        INTERVAL_VECTOR b, INTERVAL_VECTOR p);

// Krawczyk contractor
//INTERVAL_VECTOR KrawczykContractor(PIVF f, PIMF Jf, INTERVAL_VECTOR x, void* params);



#endif
