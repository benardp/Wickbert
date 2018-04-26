#ifndef VECTOR2D_H
#define VECTOR2D_H
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <sstream>

#include "scalar.h"


class Vector2d{
  friend std::ostream &operator<<(std::ostream&, const Vector2d&);
  friend std::istream &operator>>(std::istream &input, Vector2d& v);
 private:
  scalar array[2];
 public:

  Vector2d(){
    array[0] = 0;
    array[1] = 0;
  }

  Vector2d(scalar a){
    array[0] = a;
    array[1] = a;
  }

  Vector2d(scalar x, scalar y){
    array[0] = x;
    array[1] = y;
  }

  static unsigned int Dimension() { return 2; }

  /* Access Specific entries */
  scalar X() const{ return array[0]; }
  scalar Y() const{ return array[1]; }

  scalar  operator[](int i) const{ return array[i]; }
  scalar& operator[](int i)      { return array[i]; }

  /* Equality operators */
  bool operator==(const Vector2d & v) const{  return (X() == v.X()) && (Y() == v.Y()); }
  bool operator!=(const Vector2d & v) const{  return (X() != v.X()) || (Y() != v.Y()); }

  /* Dot Product */
  scalar operator*(const Vector2d & v) const {  return array[0] * v.X() + array[1] * v.Y(); }
  
  /* Vector Addition */
  Vector2d operator+(const Vector2d & v) const{  return Vector2d(array[0] + v.X(),array[1] + v.Y());}
  void operator+=(const Vector2d & v) {   array[0] += v.X();  array[1] += v.Y(); }

  /* Vector Subtraction */
  Vector2d operator-(const Vector2d & v) const{  return Vector2d( X() - v.X(), Y() - v.Y());}
  void operator-=(const Vector2d & v) { array[0] -= v.X(); array[1] -= v.Y(); }

  /* Scalar Multiplication */
  Vector2d operator*(const scalar f) const { return Vector2d(f*X(),f*Y()); }
  void operator*=(const scalar f) { array[0] *= f; array[1] *= f; }

  /* Scalar Division */
  Vector2d operator/(const scalar f) const{
    scalar f_inv = ((scalar) 1.0)/f;  
    return Vector2d(f_inv * X(), f_inv * Y());
  }

  void operator/=(const scalar f){
    scalar f_inv = ((scalar)1.0)/f;  
    array[0] *= f_inv;
    array[1] *= f_inv;
  }

  /* Euclidean Norm of this Vector */
  scalar norm() const{ return SQRT(X() * X() + Y() * Y()); }

  /* Square of the Norm of this Vector */
  scalar norm2() const{ return (X() * X() + Y() * Y()); }
   
  /* Return a normalized copy of this Vector */
  Vector2d normalized() const{ return (*this) / norm(); }
 
  scalar * loadFromArray(scalar * a){
    memcpy(array,a,Dimension() * sizeof(scalar));
    return &a[Dimension()];
  }

  scalar * writeToArray(scalar * a){
    memcpy(a,array,Dimension() * sizeof(scalar));
    return &a[Dimension()];
  }


  std::string toString() const{
    std::ostringstream s;
    s << "(" << X() << "," << Y() << ")";
    return s.str();
  } 
};


inline void normalize(Vector2d& v){
  v /= v.norm();
}

inline Vector2d operator-(const Vector2d &u){
  return Vector2d(-u.X(), -u.Y());
}


inline Vector2d operator*(scalar f, const Vector2d &v){
  return Vector2d(f * v.X(), f * v.Y());
}
  
inline std::ostream &operator<<(std::ostream &output, const Vector2d& v){
  return output << v.X() << " " << v.Y();
}

inline std::istream &operator>>(std::istream &input, Vector2d& v){
  input >> v[0];
  input >> v[1];

  return input;
}

#endif
