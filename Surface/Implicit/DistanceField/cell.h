#ifndef CELL_H
#define CELL_H

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <sstream>

#include "defs.h"
#include "vector3d.h"
#include "mdsimhash.h"

class Cell{
  friend std::ostream &operator<<(std::ostream&, const Cell&);
  friend std::istream &operator>>(std::istream &input, Cell& c);
 private:
  int location[3];
 public:
  Cell(int x = 0, int y = 0, int z = 0){
    location[0] = x;
    location[1] = y;
    location[2] = z;
  }

  Cell(const int x[3]){
    location[0] = x[0];
    location[1] = x[1];
    location[2] = x[2];
  }

  /* Access Specific entries */
  int  operator[](int i) const{
    return location[i];
  }
    
  int& operator[](int i){
    return location[i];
  }

  int X() const{
    return location[0];
  }
  
  int Y() const{
    return location[1];
  }
 
  int Z() const{
    return location[2];
  }

  Vector3d toVector3d() const{ 
    return Vector3d((scalar) X(),(scalar) Y(), (scalar) Z());
  }

  Cell operator+(const Cell &c) const { return Cell(X() + c.X(), Y() + c.Y(), Z() + c.Z()); }
  Cell operator-(const Cell &c) const { return Cell(X() - c.X(), Y() - c.Y(), Z() - c.Z()); }
  bool operator==(const Cell &c) const { return (X() == c.X()) && (Y() == c.Y()) && (Z() == c.Z());}
  bool operator!=(const Cell &c) const { return (X() != c.X()) || (Y() != c.Y()) || (Z() != c.Z());}
};  






#ifdef WIN32
struct cell_hash_compare
{
   enum
   {
      bucket_size = 4,
      min_buckets = 8
   };

   size_t operator()(const Cell& p) const
   {
      return  150001*p.X() + 157057*p.Y() + 170003*p.Z();
   }

   bool operator()(const Cell& p1, const Cell& p2) const
   {
	   if (p1.X() < p2.X())
		   return true;
	   if (p1.X() > p2.X())
		   return false;

	   if (p1.Y() < p2.Y())
		   return true;
	   if (p1.Y() > p2.Y())
		   return false;

       if (p1.Z() < p2.Z())
		   return true;   
	
	   return false;

   }
};


typedef std::set<Cell,cell_hash_compare> CellSet;
#else
struct hash_cell{
  std::hash<int> H;
  size_t operator()(const Cell& c) const{
    int bits1, bits2, bits3;
    
    bits1  = c.X() ;
    bits2  = c.Y() << 10;
    bits3  = c.Z() << 20;

    return (bits1 ^ bits2 ^ bits3);
  }
};

struct equal_cell{
  bool operator()(const Cell& a, const Cell& b) const{
    return a == b;
  }
};

typedef HASH_VERSION::hash_set< Cell, hash_cell, equal_cell > CellSet;
#endif


inline Cell toCell(const Vector3d& v, scalar cell_size){
  scalar cell_size_inverse = (scalar)1.0 / cell_size;
  return Cell( FLOOR((scalar) v.X() * cell_size_inverse),
	       FLOOR((scalar) v.Y() * cell_size_inverse),
	       FLOOR((scalar) v.Z() * cell_size_inverse));

}

inline Cell toCellInverse(const Vector3d& v, scalar cell_size_inverse){
  return Cell( FLOOR((scalar) v.X() * cell_size_inverse),
	       FLOOR((scalar) v.Y() * cell_size_inverse),
	       FLOOR((scalar) v.Z() * cell_size_inverse));
}

inline Cell toNearestCell(const Vector3d& v, scalar cell_size){
  scalar cell_size_inverse = (scalar)1.0 / cell_size;
  return Cell( ROUND((scalar) v.X() * cell_size_inverse),
	       ROUND((scalar) v.Y() * cell_size_inverse),
	       ROUND((scalar) v.Z() * cell_size_inverse));

}


inline std::ostream &operator<<(std::ostream &output, const Cell& cell){
  return output << cell.X() << " " << cell.Y() << " " << cell.Z();
}

inline std::istream &operator>>(std::istream &input, Cell& c){
  input >> c[0];
  input >> c[1];
  input >> c[2];

  return input;
}




#endif
