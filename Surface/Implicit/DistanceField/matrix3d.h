#ifndef MATRIX3D_H
#define MATRIX3D_H
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <sstream>

#include "defs.h"
#include "vector3d.h"

#define DIMENSION 3

class Matrix3d{
  friend std::ostream &operator<<(std::ostream&, const Matrix3d&);
  friend std::istream &operator>>(std::istream &input, Matrix3d& v);
 private:
  Vector3d row[3];
 public:
  /* Constructors */
  Matrix3d(){}
  Matrix3d(const Matrix3d * ptr) { memcpy(this, ptr, sizeof(Matrix3d)); }
  Matrix3d(scalar x1y1, scalar x1y2, scalar x1y3,
	   scalar x2y1, scalar x2y2, scalar x2y3,
	   scalar x3y1, scalar x3y2, scalar x3y3);  
  Matrix3d(const Vector3d& r0, const Vector3d& r1, const Vector3d& r2);

  static Matrix3d I();
  static Matrix3d Zero();

  std::string toString() const;
  
  /* Matrix Vector Product */
  Vector3d operator*(const Vector3d & v) const;
  
  /* Matrix Addition */
  Matrix3d  operator+ (const Matrix3d & m) const;
  Matrix3d& operator+=(const Matrix3d & m);
  
  /* Matrix Multiplication */
  Matrix3d operator*(const Matrix3d & m) const;

  /* Matrix Subtraction */
  Matrix3d  operator- (const Matrix3d & m) const;
  Matrix3d& operator-=(const Matrix3d & m);

  /* Scalar Multiplication */
  Matrix3d  operator* (scalar f) const;
  Matrix3d& operator*=(scalar f);

  /* Scalar Division */
  Matrix3d  operator/(scalar f) const;
  Matrix3d& operator/=(scalar f);

  /* The Trace of this Matrix */
  scalar trace() const;

  Matrix3d transpose() const;
  scalar det() const;

  Matrix3d adjoint() const;
  Matrix3d inverse() const;

  /* Accessors */
  Vector3d&       operator[](int i)       { return row[i]; }
  const Vector3d& operator[](int i) const { return row[i]; }

  scalar& operator()(int i, int j)       { return row[i][j]; }
  scalar  operator()(int i, int j) const { return row[i][j]; }

  Vector3d col(int i) const {return Vector3d(row[0][i],row[1][i],row[2][i]);} 
};


/*
 * Non-member prototypes
 */
void transpose(Matrix3d & m);
Matrix3d outer(const Vector3d &u, const Vector3d &v);
Matrix3d outer(const Vector3d &u);
Matrix3d operator*(scalar f, const Matrix3d &m);
Matrix3d Star(const Vector3d &a);



/* 
 * Constructors
 */
inline Matrix3d::Matrix3d(scalar x1y1, scalar x1y2, scalar x1y3,
			  scalar x2y1, scalar x2y2, scalar x2y3,
			  scalar x3y1, scalar x3y2, scalar x3y3){
  row[0][0] = x1y1;
  row[0][1] = x1y2;
  row[0][2] = x1y3;

  row[1][0] = x2y1;
  row[1][1] = x2y2;
  row[1][2] = x2y3;

  row[2][0] = x3y1;
  row[2][1] = x3y2;
  row[2][2] = x3y3;
}


inline Matrix3d::Matrix3d(const Vector3d& r0, const Vector3d& r1, const Vector3d& r2){
  row[0] = r0;
  row[1] = r1;
  row[2] = r2;
}


/*
 * Identity matrices
 */
inline Matrix3d Matrix3d::I(){
  return Matrix3d(1.0, 0.0, 0.0,
		  0.0, 1.0, 0.0,
		  0.0, 0.0, 1.0);
}

inline Matrix3d Matrix3d::Zero(){
  return Matrix3d(0.0, 0.0, 0.0,
		  0.0, 0.0, 0.0,
		  0.0, 0.0, 0.0);
}



/* A simple text description */      
inline std::string Matrix3d::toString() const{
  std::ostringstream s;
  for(int i = 0; i < DIMENSION; i++){
    s << row[i] << std::endl;
  }
  return s.str();
}
  


/* Matrix Vector Product */
inline Vector3d Matrix3d::operator*(const Vector3d & v) const{
  return Vector3d(row[0] * v, row[1] * v, row[2] * v);
}

/* Matrix Addition */
inline Matrix3d Matrix3d::operator+(const Matrix3d & m) const{
  Matrix3d result(this);

  result[0] += m[0];
  result[1] += m[1];
  result[2] += m[2];

  return result;
}

inline Matrix3d& Matrix3d::operator+=(const Matrix3d & m){
  row[0] += m[0];
  row[1] += m[1];
  row[2] += m[2];

  return *this;
}


/* Matrix Subtraction */
inline Matrix3d Matrix3d::operator-(const Matrix3d & m) const{
  Matrix3d result(this);

  result[0] -= m[0];
  result[1] -= m[1];
  result[2] -= m[2];

  return result;
}

inline Matrix3d& Matrix3d::operator-=(const Matrix3d & m){
  row[0] -= m[0];
  row[1] -= m[1];
  row[2] -= m[2];

  return *this;
}


/* Scalar Multiplication */
inline Matrix3d Matrix3d::operator*(const scalar f) const{
  Matrix3d result(this);

  result[0] *= f;
  result[1] *= f;
  result[2] *= f;

  return result;
}

inline Matrix3d& Matrix3d::operator*=(const scalar f){
  row[0] *= f;
  row[1] *= f;
  row[2] *= f;

  return *this;
}

/* Scalar Division */
inline Matrix3d Matrix3d::operator/(const scalar f) const{
  Matrix3d result(this);
  scalar f_inv = (scalar) 1.0 / f;

  result[0] *= f_inv;
  result[1] *= f_inv;
  result[2] *= f_inv;

  return result;
}

inline Matrix3d& Matrix3d::operator/=(scalar f){
  Matrix3d result(this);
  scalar f_inv = (scalar) 1.0 / f;

  row[0] *= f_inv;
  row[1] *= f_inv;
  row[2] *= f_inv;

  return *this;
}


/* Matrix Multiplication */
inline Matrix3d Matrix3d::operator*(const Matrix3d &m) const{
  const Matrix3d mt = m.transpose();
 
  return Matrix3d(row[0] * mt[0], row[0] * mt[1], row[0] * mt[2],
		  row[1] * mt[0], row[1] * mt[1], row[1] * mt[2],
		  row[2] * mt[0], row[2] * mt[1], row[2] * mt[2]);
}


/* Matrix determinate */
inline scalar Matrix3d::det() const{
  return row[0] * cross(row[1],row[2]);
}


/* trace */
inline scalar Matrix3d::trace() const{
  return row[0][0] + row[1][1] + row[2][2];
}


/* Adjoint and Inverse */
inline Matrix3d Matrix3d::adjoint() const{
  return Matrix3d(cross(row[1],row[2]),
		  cross(row[2],row[0]),
		  cross(row[0],row[1]));
}


inline Matrix3d Matrix3d::inverse() const{
  Matrix3d adj(adjoint());

  scalar det = adj[0] * row[0];

  adj = adj.transpose();
  adj /= det;

  return adj;
}

/* transpose */
inline Matrix3d Matrix3d::transpose() const{
  return Matrix3d(row[0][0],row[1][0],row[2][0],
		  row[0][1],row[1][1],row[2][1],
		  row[0][2],row[1][2],row[2][2]);
}



/* Non member functions */
inline void transpose2(Matrix3d & m){
  scalar temp;

  temp = m[0][1];
  m[0][1] = m[1][0];
  m[1][0] = temp;

  temp = m[0][2];
  m[0][2] = m[2][0];
  m[2][0] = temp;

  temp = m[1][2];
  m[1][2] = m[2][1];
  m[2][1] = temp;
}


	
inline Matrix3d outer(const Vector3d &u, const Vector3d &v){
  return Matrix3d(u.X() * v.X(), u.X() * v.Y(), u.X() * v.Z(),
		  u.Y() * v.X(), u.Y() * v.Y(), u.Y() * v.Z(),
		  u.Z() * v.X(), u.Z() * v.Y(), u.Z() * v.Z());
}

inline Matrix3d outer(const Vector3d &u){
  return Matrix3d(u.X() * u.X(), u.X() * u.Y(), u.X() * u.Z(),
		  u.Y() * u.X(), u.Y() * u.Y(), u.Y() * u.Z(),
		  u.Z() * u.X(), u.Z() * u.Y(), u.Z() * u.Z());
}


inline Matrix3d operator*(scalar f, const Matrix3d &m){
  Matrix3d p = m;

  p *= f;

  return p;
}
  


inline std::ostream &operator<<(std::ostream &output, const Matrix3d& m){
  output << m[0] << " ";
  output << m[1] << " ";
  output << m[2];

  return output;
}

inline std::istream &operator>>(std::istream &input, Matrix3d& m){
  input >> m[0];
  input >> m[1];
  input >> m[2];

  return input;
}

#endif
