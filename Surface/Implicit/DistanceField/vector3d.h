#ifndef VECTOR3D_H
#define VECTOR3D_H
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <sstream>

#include "scalar.h"


class Vector3d{
  friend std::ostream &operator<<(std::ostream&, const Vector3d&);
  friend std::istream &operator>>(std::istream &input, Vector3d& v);
 private:
  scalar array[3];
 public:

  Vector3d(){
    array[0] = 0;
    array[1] = 0;
    array[2] = 0;
  }

  Vector3d(scalar a){
    array[0] = a;
    array[1] = a;
    array[2] = a;
  }

  Vector3d(scalar x, scalar y, scalar z){
    array[0] = x;
    array[1] = y;
    array[2] = z;
  }

  static unsigned int Dimension() { return 3; }

  /* Access Specific entries */
  scalar X() const{ return array[0]; }
  scalar Y() const{ return array[1]; }
  scalar Z() const{ return array[2]; }

  scalar  operator[](int i) const{ return array[i]; }
  scalar& operator[](int i)      { return array[i]; }

  /* Equality operators */
  bool operator==(const Vector3d & v) const{  return (X() == v.X()) && (Y() == v.Y()) && (Z() == v.Z()); }
  bool operator!=(const Vector3d & v) const{  return (X() != v.X()) || (Y() != v.Y()) || (Z() != v.Z()); }

  /* Dot Product */
  scalar operator*(const Vector3d & v) const {  return array[0] * v.X() + array[1] * v.Y() + array[2] * v.Z(); }
  
  /* Vector Addition */
  Vector3d operator+(const Vector3d & v) const{  return Vector3d(array[0] + v.X(),array[1] + v.Y(),array[2] + v.Z());}
  void operator+=(const Vector3d & v) {   array[0] += v.X();  array[1] += v.Y();  array[2] += v.Z(); }

  /* Vector Subtraction */
  Vector3d operator-(const Vector3d & v) const{  return Vector3d( X() - v.X(), Y() - v.Y(), Z() - v.Z());}
  void operator-=(const Vector3d & v) { array[0] -= v.X(); array[1] -= v.Y(); array[2] -= v.Z(); }

  /* Scalar Multiplication */
  Vector3d operator*(const scalar f) const { return Vector3d(f*X(),f*Y(),f*Z()); }
  void operator*=(const scalar f) { array[0] *= f; array[1] *= f; array[2] *= f; }

  /* Scalar Division */
  Vector3d operator/(const scalar f) const{
    scalar f_inv = ((scalar) 1.0)/f;  
    return Vector3d(f_inv * X(), f_inv * Y(), f_inv * Z());
  }

  void operator/=(const scalar f){
    scalar f_inv = ((scalar)1.0)/f;  
    array[0] *= f_inv;
    array[1] *= f_inv;
    array[2] *= f_inv;
  }

  /* Euclidean Norm of this Vector */
  scalar norm() const{ return SQRT(X() * X() + Y() * Y() + Z() * Z()); }

  /* Square of the Norm of this Vector */
  scalar norm2() const{ return (X() * X() + Y() * Y() + Z() * Z()); }
   
  /* Return a normalized copy of this Vector */
  Vector3d normalized() const{ return (*this) / norm(); }
 
  /* Cross Product */
  Vector3d cross(const Vector3d & v) const{
    return Vector3d(    Y() * v.Z() - v.Y() *   Z(),
		      v.X() *   Z() -   X() * v.Z(),
		        X() * v.Y() - v.X() *   Y());
  }


  scalar * loadFromArray(scalar * a){
    memcpy(array,a,3 * sizeof(scalar));
    return &a[3];
  }

  scalar * writeToArray(scalar * a){
    memcpy(a,array,3 * sizeof(scalar));
    return &a[3];
  }


  std::string toString() const{
    std::ostringstream s;
    s << "(" << X() << "," << Y() << "," << Z() << ")";
    return s.str();
  } 

};


inline void normalize(Vector3d& v){
  v /= v.norm();
}

inline Vector3d operator-(const Vector3d &u){
  return Vector3d(-u[0], -u[1], -u[2]);
}

/* Cross Product */
inline Vector3d cross(const Vector3d & u, const Vector3d &v){
  return Vector3d( u[1] * v[2] - v[1] * u[2],
		   v[0] * u[2] - u[0] * v[2],
		   u[0] * v[1] - v[0] * u[1] );
}


inline Vector3d operator*(scalar f, const Vector3d &v){
  return Vector3d(f * v.X(), f * v.Y(), f * v.Z());
}
  
inline std::ostream &operator<<(std::ostream &output, const Vector3d& v){
  return output << v.X() << " " << v.Y() << " " << v.Z();
}

inline std::istream &operator>>(std::istream &input, Vector3d& v){
  input >> v[0];
  input >> v[1];
  input >> v[2];

  return input;
}

#endif
