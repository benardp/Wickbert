#ifndef PARTICLE_H
#define PARTICLE_H
 
#include "vector3d.h"

class Particle{
  friend std::ostream &operator<<(std::ostream &output, const Particle&);
  friend std::istream &operator>>(std::istream &input,        Particle& p);
 public:
  Vector3d location;
  scalar   radius;

 public: 
  Particle() : location(Vector3d(0,0,0)), radius(1.0) {}
  Particle(const Vector3d& loc){ location = loc; radius = 1.0;}
  Particle(const Vector3d& loc, const scalar rad){ location = loc; radius = rad; }
   
  const Vector3d&  Location() const{ return location; }
  void  SetLocation(const Vector3d& new_location) { location = new_location; }

  const scalar&  Radius() const{ return radius; }
  void  SetRadius(const scalar& new_radius) { radius = new_radius; }
};



inline std::ostream &operator<<(std::ostream &output, const Particle& p){
  output << p.location << " ";
  output << p.radius;
  return output;
}

inline std::istream &operator>>(std::istream &input, Particle& p){
  input >> p.location; 
  input >> p.radius; 

  return input;
}


#endif
