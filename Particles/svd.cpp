/**
 * Implementation of singular value decomposition.
 * @file svd.cpp
 * @date July 9, 2001
 * @author Ed Bachta
 * @remarks Based on code from Numerical Recipies
 */


#include "svd.h"


/** 
 * Computes sqrt(a^2 + b^2 ). No destructive underflow or overflow.
 */
double SVD::pythag(double a, double b) { 

  double absa,absb; 

  absa=fabs(a); 
  absb=fabs(b); 

  if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa)); 
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb))); 

} // end pythag


/**
 * Fixes ill-conditioning by removing small singular values.
 * @param w Singular value vector to be conditioned.
 */
void SVD::condition(TNT::Vector<double> &w) {

  int i, n = w.size();
  double max_v = 0.0;

  TNT::Vector<double> w_mag(n);

  // Compute magnitudes and find largest
  for(i=1;i<=n;i++) {
    w_mag(i) = fabs(w(i));
    if (w_mag(i) > max_v) max_v = w_mag(i);    
  }

  if (max_v==0.0) return; // No non-zero singular values

  // Find smallest magnitude
  int min_i = 0;
  double min_v = max_v;
  for(i=1;i<=n;i++) {
    if ((w_mag(i) != 0.0) && (w_mag(i) <= min_v)) {
      min_v = w_mag(i);
      min_i = i;
    }
  }

  //if (min_i==0) return; // No non-zero singular values

  // Test condition number
  while ((min_v / max_v) < 0.000001) { 

    std::cerr << "SVD Condition number: " << min_v/max_v << std::endl;

    // Zero out smallest singular value
    w(min_i) = w_mag(min_i) = 0.0;

    // Find next smallest singular value
    min_i = 0;
    min_v = max_v;
    for(i=1;i<=n;i++) {
      if ((w_mag(i) != 0.0) && (w_mag(i) < min_v)) {
        min_v = w_mag(i);
        min_i = i;
      }
    }

    if (min_i==0) return; // No more non-zero singular values

  }

} // end condition


/**
 * Computes SVD using TNT Matrices.
 * @param a Matrix to decompose (replaced by orthogonal matrix U).
 * @param w Vector of singular values (actually a diagonal matrix).
 * @param v Orthogonal matrix V.
 * @note Orthogonal matrix U is returned in a!
 * @returns False if convergance fails.
 */
bool SVD::compute(TNT::Matrix<double> &a,
                  TNT::Vector<double> &w,
                  TNT::Matrix<double> &v) {

  int m = a.num_rows();
  int n = a.num_cols();

  int flag,i,its,j,jj,k,l,nm; 
  double anorm,c,f,g,h,s,scale,x,y,z;

  TNT::Vector<double> rv1(n);

  g=scale=anorm=0.0; 
 
  for (i=1;i<=n;i++) { 
    l=i+1; 
    rv1(i)=scale*g; 
    g=s=scale=0.0; 
    if (i <= m) { 
      for (k=i;k<=m;k++) scale += fabs(a(k,i)); 
      if (scale) { 
        for (k=i;k<=m;k++) { 
          a(k,i) /= scale; 
          s += a(k,i)*a(k,i); 
        } 
        f=a(i,i); 
        g = -SIGN(sqrt(s),f); 
        h=f*g-s; 
        a(i,i)=f-g; 
        for (j=l;j<=n;j++) { 
          for (s=0.0,k=i;k<=m;k++) s += a(k,i)*a(k,j); 
          f=s/h; 
          for (k=i;k<=m;k++) a(k,j) += f*a(k,i); 
        } 
        for (k=i;k<=m;k++) a(k,i) *= scale; 
      } 
    } 
    w(i)=scale *g; 
    g=s=scale=0.0; 
    if (i <= m && i != n) { 
      for (k=l;k<=n;k++) scale += fabs(a(i,k)); 
      if (scale) {
        for (k=l;k<=n;k++) { 
          a(i,k) /= scale; 
          s += a(i,k)*a(i,k); 
        } 
        f=a(i,l); 
        g = -SIGN(sqrt(s),f); 
        h=f*g-s; 
        a(i,l)=f-g; 
        for (k=l;k<=n;k++) rv1(k)=a(i,k)/h; 
        for (j=l;j<=m;j++) { 
          for (s=0.0,k=l;k<=n;k++) s += a(j,k)*a(i,k); 
          for (k=l;k<=n;k++) a(j,k) += s*rv1(k); 
        } 
        for (k=l;k<=n;k++) a(i,k) *= scale; 
      } 
    } 
    anorm=FMAX(anorm,(fabs(w(i))+fabs(rv1(i)))); 
  } 
  for (i=n;i>=1;i--) { 
    if (i < n) { 
      if (g) { 
        for (j=l;j<=n;j++) v(j,i)=(a(i,j)/a(i,l))/g; 
        for (j=l;j<=n;j++) { 
          for (s=0.0,k=l;k<=n;k++) s += a(i,k)*v(k,j); 
          for (k=l;k<=n;k++) v(k,j) += s*v(k,i); 
        } 
      } 
      for (j=l;j<=n;j++) v(i,j)=v(j,i)=0.0; 
    } 
    v(i,i)=1.0; 
    g=rv1(i); 
    l=i; 
  } 
  for (i=IMIN(m,n);i>=1;i--) { 
    l=i+1; 
    g=w(i); 
    for (j=l;j<=n;j++) a(i,j)=0.0; 
    if (g) { 
      g=1.0/g; 
      for (j=l;j<=n;j++) { 
        for (s=0.0,k=l;k<=m;k++) s += a(k,i)*a(k,j); 
        f=(s/a(i,i))*g; 
        for (k=i;k<=m;k++) a(k,j) += f*a(k,i); 
      } 
      for (j=i;j<=m;j++) a(j,i) *= g; 
    } else for (j=i;j<=m;j++) a(j,i)=0.0; 
    ++a(i,i); 
  } 
  for (k=n;k>=1;k--) { 
    for (its=1;its<=limit;its++) { 
      flag=1; for (l=k;l>=1;l--) { 
        nm=l-1; 
        if ((double)(fabs(rv1(l))+anorm) == anorm) { 
          flag=0; 
          break; 
        } 
        if ((double)(fabs(w(nm))+anorm) == anorm) break; 
      } 
      if (flag) { 
        c=0.0;  
        s=1.0; 
        for (i=l;i<=k;i++) {
          f=s*rv1(i); 
          rv1(i)=c*rv1(i); 
          if ((double)(fabs(f)+anorm) == anorm) break; 
          g=w(i); 
          h=pythag(f,g); 
          w(i)=h; 
          h=1.0/h; 
          c=g*h; 
          s = -f*h; 
          for (j=1;j<=m;j++) { 
            y=a(j,nm); 
            z=a(j,i); 
            a(j,nm)=y*c+z*s; 
            a(j,i)=z*c-y*s; 
          } 
        } 
      } 
      z=w(k); 
      if (l == k) { 
        if (z < 0.0) { 
          w(k) = -z; 
          for (j=1;j<=n;j++) v(j,k) = -v(j,k); 
        } 
        break; 
      } 
      if (its == limit) return false; 
      x=w(l);
      nm=k-1; 
      y=w(nm); 
      g=rv1(nm); 
      h=rv1(k); 
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y); 
      g=pythag(f,1.0); 
      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x; 
      c=s=1.0; 
      for (j=l;j<=nm;j++) { 
        i=j+1; 
        g=rv1(i); 
        y=w(i); 
        h=s*g;
        g=c*g; 
        z=pythag(f,h); 
        rv1(j)=z; 
        c=f/z; 
        s=h/z; 
        f=x*c+g*s; 
        g = g*c-x*s; 
        h=y*s; 
        y *= c; 
        for (jj=1;jj<=n;jj++) { 
          x=v(jj,j); 
          z=v(jj,i); 
          v(jj,j)=x*c+z*s; 
          v(jj,i)=z*c-x*s; 
        } 
        z=pythag(f,h); 
        w(j)=z;
        if (z) { 
          z=1.0/z; 
          c=f*z; 
          s=h*z; 
        } 
        f=c*g+s*y; 
        x=c*y-s*g;
        for (jj=1;jj<=m;jj++) { 
          y=a(jj,j); 
          z=a(jj,i); 
          a(jj,j)=y*c+z*s; 
          a(jj,i)=z*c-y*s; 
        } 
      } 
      rv1(l)=0.0; 
      rv1(k)=f; 
      w(k)=x; 
    } 
  } 

  return true;

} // end compute


/**
 * Solves a system of equations using SVD components.
 * @param u Orthogonal SVD matrix U.
 * @param w Singular value vector W.
 * @param v Orthogonal SVD matrix V.
 * @param b Solution vector.
 * @param x Variable vector.
 */
void SVD::solve(TNT::Matrix<double> &u,
                TNT::Vector<double> &w,
                TNT::Matrix<double> &v,
                TNT::Vector<double> b,
                TNT::Vector<double> &x) {
  
  int m = u.num_rows();
  int n = u.num_cols();

  int j,i; 
  double s;

  TNT::Vector<double> tmp(n); 
 
  for (j=1;j<=n;j++) { 
    s=0.0; 
    if (w(j)) { 
      // Nonzero result only if wj is nonzero. 
      for (i=1;i<=m;i++) s += u(i,j)*b(i); 
      s /= w(j); // This is the divide by wj . 
    } 
    tmp(j)=s; 
  } 

  x = v*tmp;

} // end solve


/**
 * Solves a system of equations.
 * @param A Matrix of coefficients.
 * @param x Vector of variables.
 * @param b Solution vector.
 * @returns False if SVD did not converge, true otherwise
 */
bool SVD::solve(TNT::Matrix<double> A,
                TNT::Vector<double> &x,
                TNT::Vector<double> b) {

  int m = A.num_rows();
  int n = A.num_cols();

  TNT::Matrix<double> U(m,n);
  TNT::Matrix<double> V(n,n);
  TNT::Vector<double> W(n);
  
  U = A;

  if (!compute(U, W, V)) return false;

  condition(W);

  solve(U, W, V, b, x);

  return true;

} // end solve



