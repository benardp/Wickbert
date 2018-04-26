#include "scalarlattice.h"

const Cell ScalarLattice::N6i[6] = {
  Cell(-1, 0, 0),
  Cell( 1, 0, 0),
  Cell( 0,-1, 0),
  Cell( 0, 1, 0),
  Cell( 0, 0,-1),
  Cell( 0, 0, 1)
};



/*
 * Determine the normal at this point using centered differences
 */
Vector3d ScalarLattice::Normal(const Cell& c) const{
  scalar c_val = Value(c);

  scalar dx;
  scalar dy;
  scalar dz;

  scalar forward, backward;

  forward = Value(c + Cell(1,0,0)) - c_val;
  backward = c_val - Value(c - Cell(1,0,0));
  
  dx = (forward + backward) / 2;
  
  forward = Value(c + Cell(0,1,0)) - c_val;
  backward = c_val - Value(c - Cell(0,1,0));

  dy = (forward + backward) / 2;

  forward = Value(c + Cell(0,0,1)) - c_val;
  backward = c_val - Value(c - Cell(0,0,1));

  dz = (forward + backward) / 2;

  return InverseVoxelSize() * Vector3d(dx,dy,dz);
}


/*
 * Using the supplied vector V, find the product V*N using upwind differencing
 */
scalar ScalarLattice::Upwind(const Cell& c, const Vector3d& V) const{
  return UpwindNormal(c,V) * V;
}


/*
 * Using the supplied vector V, find the Normal using upwind differencing
 *
 * Assumes all 6 neighbors of c have valid values
 */
Vector3d ScalarLattice::UpwindNormal(const Cell& c, const Vector3d& V) const{
  scalar c_val = Value(c);

  scalar dx;
  scalar dy;
  scalar dz;

  if(V.X() > 0.0){
    dx = c_val - Value(c - Cell(1,0,0));
  } else {
    dx = Value(c + Cell(1,0,0)) - c_val;
  }

  if(V.Y() > 0.0){
    dy = c_val - Value(c - Cell(0,1,0));
  } else {
    dy = Value(c + Cell(0,1,0)) - c_val;
  }

  if(V.Z() > 0.0){
    dz = c_val - Value(c - Cell(0,0,1));
  } else {
    dz = Value(c + Cell(0,0,1)) - c_val;
  }

  return InverseVoxelSize() * Vector3d(dx,dy,dz);
}


#define phi(i,j,k)  p[1+i][1+j][1+k]
//#define SLOW_METHOD
scalar ScalarLattice::MeanCurvature(const Cell& c) const{
  scalar mcurve  = 0.0;

#ifdef SLOW_METHOD
  scalar p[3][3][3];

  for(int i = -1; i <= 1; i++){
    for(int j = -1; j <= 1; j++){
      for(int k = -1; k <= 1; k++){
	phi(i,j,k) = Value(c + Cell(i,j,k));
      }
    }
  }

  scalar inverse = InverseVoxelSize();
  scalar inverse_squared = inverse*inverse;
  
  //Centered difference approx. to 1st derivatives
  scalar x = 0.5*inverse * (phi(1,0,0) - phi(-1,0,0));
  scalar y = 0.5*inverse * (phi(0,1,0) - phi(0,-1,0));
  scalar z = 0.5*inverse * (phi(0,0,1) - phi(0,0,-1));

  //Squares of 1st derivatives
  scalar x2 = x*x;
  scalar y2 = y*y;
  scalar z2 = z*z;

  //Centered difference approx. to 2nd derivatives
  scalar xx = inverse_squared * (phi(1,0,0) - 2*phi(0,0,0) + phi(-1,0,0));
  scalar yy = inverse_squared * (phi(0,1,0) - 2*phi(0,0,0) + phi(0,-1,0));
  scalar zz = inverse_squared * (phi(0,0,1) - 2*phi(0,0,0) + phi(0,0,-1));

  //Mixed partials
  scalar xy = 0.25*inverse * (phi(1,1,0) - phi(1,-1,0) - phi(-1,1,0) + phi(-1,-1,0));
  scalar xz = 0.25*inverse * (phi(1,0,1) - phi(1,0,-1) - phi(-1,0,1) + phi(-1,0,-1));
  scalar yz = 0.25*inverse * (phi(0,1,1) - phi(0,1,-1) - phi(0,-1,1) + phi(0,-1,-1));
  
  mcurve = (yy + zz)*x2 + (xx + zz)*y2 + (xx + yy)*z2;
  mcurve -= 2*(x*y*xy + x*z*xz + y*z*yz);

  mcurve /= (x2 + y2 + z2)*SQRT(x2 + y2 + z2);
#else
  for(int i = 0; i < 6; i++){
    mcurve += Value(c + N6i[i]);
  }  
  mcurve -= 6*Value(c);  
  mcurve *= InverseVoxelSize() * InverseVoxelSize();

#endif
  return mcurve;
}


Cell ScalarLattice::ToCell(const Vector3d& v) const{
  Vector3d vscaled = InverseVoxelSize() * v;

  return Cell((int)FLOOR(vscaled.X()),
	      (int)FLOOR(vscaled.Y()),
	      (int)FLOOR(vscaled.Z()));  
}

scalar ScalarLattice::InterpValue(const Vector3d& v) const{
  Vector3d vscaled = InverseVoxelSize() * v;

  Cell base((int)FLOOR(vscaled.X()),
	    (int)FLOOR(vscaled.Y()),
	    (int)FLOOR(vscaled.Z()));

  scalar r = vscaled.X() - base.X();
  scalar s = vscaled.Y() - base.Y();
  scalar t = vscaled.Z() - base.Z();

  scalar result = 0.0;

  for(int x = 0; x <= 1; x++){
    scalar coef_r = ((x == 0) ? (1-r) : r);
    for(int y = 0; y <= 1; y++){
      scalar coef_s = ((y == 0) ? (1-s) : s);
      for(int z = 0; z <= 1; z++){	
	scalar coef_t = ((z == 0) ? (1-t) : t);
	scalar coef = coef_r * coef_s * coef_t;
	result += coef * Value(base + Cell(x,y,z));
      }
    }
  }

  return result;
}


Vector3d ScalarLattice::InterpNormal(const Vector3d& v) const{
  Vector3d vscaled = InverseVoxelSize() * v;

  Cell base((int)FLOOR(vscaled.X()),
	    (int)FLOOR(vscaled.Y()),
	    (int)FLOOR(vscaled.Z()));
  
  scalar r = vscaled.X() - base.X();
  scalar s = vscaled.Y() - base.Y();
  scalar t = vscaled.Z() - base.Z();

  Vector3d result(0,0,0);

  for(int x = 0; x <= 1; x++){
    scalar coef_r = ((x == 0) ? (1-r) : r);
    for(int y = 0; y <= 1; y++){
      scalar coef_s = ((y == 0) ? (1-s) : s);
      for(int z = 0; z <= 1; z++){	
	scalar coef_t = ((z == 0) ? (1-t) : t);
	scalar coef = coef_r * coef_s * coef_t;
	result += coef * Normal(base + Cell(x,y,z));
      }
    }
  }

  return result;
}

scalar ScalarLattice::InterpMeanCurvature(const Vector3d& v) const{
  Vector3d vscaled = InverseVoxelSize() * v;

  Cell base((int)FLOOR(vscaled.X()),
	    (int)FLOOR(vscaled.Y()),
	    (int)FLOOR(vscaled.Z()));

  scalar r = vscaled.X() - base.X();
  scalar s = vscaled.Y() - base.Y();
  scalar t = vscaled.Z() - base.Z();

  scalar result = 0.0;

  for(int x = 0; x <= 1; x++){
    scalar coef_r = ((x == 0) ? (1-r) : r);
    for(int y = 0; y <= 1; y++){
      scalar coef_s = ((y == 0) ? (1-s) : s);
      for(int z = 0; z <= 1; z++){	
	scalar coef_t = ((z == 0) ? (1-t) : t);
	scalar coef = coef_r * coef_s * coef_t;
	result += coef * MeanCurvature(base + Cell(x,y,z));
      }
    }
  }

  return result;
}


