#include "utils.h"
#include "triangle.h"

/*
 * Ratio of the inradius to the circumradius (scaled by 3) of a tetrahedron
 *
 * Formula from:
 *  "Relationship between Tetrahedron Shape Measures"
 *   A. Liu and B. Joe
 *   BIT vol. 34(1994) 268--287
 *
 *
 *  Returns value in the range [0,1] which corresponds to degenerate 
 *  and regular tetrahedra respectively.
 * 
 */
scalar TetrahedronRadiusRatio(const Vector3d& v0, const Vector3d& v1, const Vector3d& v2, const Vector3d& v3){
  //The signed tet volume
  const scalar volume = ((v1 - v0) * cross(v2 - v0, v3 - v0)) / 6.0;

  //products of opposite edges
  const scalar a = (v1 - v0).norm() * (v3 - v2).norm(); // 01-23
  const scalar b = (v2 - v0).norm() * (v3 - v1).norm(); // 02-13
  const scalar c = (v3 - v0).norm() * (v2 - v1).norm(); // 03-12

  //Sum of the areas of the 4 faces
  const scalar surface_area = 
    Triangle(v0,v1,v2).SurfaceArea() + 
    Triangle(v0,v2,v3).SurfaceArea() + 
    Triangle(v0,v1,v3).SurfaceArea() + 
    Triangle(v1,v2,v3).SurfaceArea(); 

  const scalar numerator = 216 * volume * volume;
  const scalar denominator = surface_area * SQRT((a + b + c)*(a + b - c)*(a + c - b)*(b + c - a));

  return numerator/denominator;
}



std::string GetExtension(const std::string& str){
  std::string::size_type pos = str.rfind(".");
  
  if(pos == std::string::npos){ return "";}
  else { return str.substr(pos+1,str.length() - pos); }
}

std::string GetNonExtension(const std::string& str){
  std::string::size_type pos = str.rfind(".");
  
  if(pos == std::string::npos){ return str;}
  else { return str.substr(0,pos); }
}




bool TestUtils(){
  Vector3d regular_tet[4] = { Vector3d( SQRT(3.0)/3.0,    0, 0),
                              Vector3d(-SQRT(3.0)/6.0,  0.5, 0),
                              Vector3d(-SQRT(3.0)/6.0, -0.5, 0),
                              Vector3d(             0,    0, SQRT(6.0)/3.0) };
  Vector3d degenerate_tet[4] = { Vector3d(  0,   0, 0),
				 Vector3d(  1,   0, 0),
				 Vector3d(  0,   1, 0),
				 Vector3d(0.5, 0.5, 0) };

  std::cout << "utils.h  TetrahedronRadiusRatio():";
  scalar reg_value = TetrahedronRadiusRatio(regular_tet[0],regular_tet[1],regular_tet[2],regular_tet[3]);
  scalar deg_value = TetrahedronRadiusRatio(degenerate_tet[0],degenerate_tet[1],degenerate_tet[2],degenerate_tet[3]);


  assert(FABS(reg_value - 1.0) < 1e-4);
  assert(FABS(deg_value - 0.0) < 1e-4);

  std::cout << " passed" << std::endl;  


  std::cout << "utils.h  GetExtension():";
  assert(GetExtension("hello.1") == "1");
  assert(GetExtension("hello11") == "");
  assert(GetExtension("hello.1.2") == "2");
  assert(GetExtension("hello.") == "");
  assert(GetExtension(".") == "");
  assert(GetExtension("1.") == "");
  assert(GetExtension("") == "");
  std::cout << " passed" << std::endl;  

  std::cout << "utils.h  GetNonExtension():";
  assert(GetNonExtension("hello.1") == "hello");
  assert(GetNonExtension("hello11") == "hello11");
  assert(GetNonExtension("hello.1.2") == "hello.1");
  assert(GetNonExtension("hello.") == "hello");
  assert(GetNonExtension(".") == "");
  assert(GetNonExtension("1.") == "1");
  assert(GetNonExtension("") == "");
  std::cout << " passed" << std::endl;  

  return true;
}
