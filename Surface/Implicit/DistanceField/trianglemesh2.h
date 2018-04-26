#ifndef TRIANGLEMESH2_H
#define TRIANGLEMESH2_H

#include <vector>

//Are all of these necessary? I don't see the point... - Elmar
//#include "renderable.h"
//#include "sparsescalarlattice.h"
//#include "boundingvolume.h"

#include "triangle.h"
#include "index.h"
//#include "matrix3d.h" - just referenced not used - Elmar
class Matrix3d;



class IndexTriple{
  friend std::ostream &operator<<(std::ostream &output, const IndexTriple&);
  friend std::istream &operator>>(std::istream &input,        IndexTriple&);
 private:
  Index triple[3];
 public:
  IndexTriple() { triple[0] = 0; triple[1] = 0; triple[2] = 0;}
  IndexTriple(Index a,Index b,Index c) { triple[0] = a; triple[1] = b; triple[2] = c;}
  
  const Index& A() const {return triple[0];}
        Index& A()       {return triple[0];}
  const Index& B() const {return triple[1];}
        Index& B()       {return triple[1];}
  const Index& C() const {return triple[2];}
        Index& C()       {return triple[2];}

  const Index& operator[](unsigned int i) const{ return triple[i]; }
        Index& operator[](unsigned int i)      { return triple[i]; }
};

inline std::ostream &operator<<(std::ostream &output, const IndexTriple& t){
  output << t.A() << " " << t.B() << " " << t.C();
  return output;
}

inline std::istream &operator>>(std::istream &input,        IndexTriple& t){
  input >> t.A() >> t.B() >> t.C();
  return input;
}



class TriangleMesh2 : public Renderable {
  friend std::ostream &operator<<(std::ostream &output, const TriangleMesh2&);
  friend std::istream &operator>>(std::istream &input,        TriangleMesh2&);
 public:
  std::vector< Vector3d > vertices;
  std::vector< IndexTriple > faces;

  std::vector< Vector3d > normals;
  std::vector< IndexTriple > face_normals;

 public:
  TriangleMesh2(){}
  TriangleMesh2(const std::string filename);
  virtual ~TriangleMesh2();

  Triangle operator[](unsigned int i) const;
  Index size() const { return faces.size(); }
  bool IsLoaded() { return loadedFlag; }

  scalar Volume() const;
  void MassProperties (scalar& volume, Vector3d& cm, Matrix3d& inertia) const;

  static TriangleMesh2 UnitCube();

  void Translate(const Vector3d& v);
  void Scale(const Vector3d& v);
  void Transform(const Matrix3d& m);


  AABoundingBox BoundingBox() const;
/*   void Intersect(const Ray& r, std::vector< scalar >& intersections) const; */

  bool LoadFromOBJFile(const std::string& filename);
  void LoadFromSMESHFile(const std::string& filename);

  void SaveToOBJFile(const std::string& filename);

  virtual void Render();

private:

  bool loadedFlag;
};




inline std::ostream &operator<<(std::ostream &output, const TriangleMesh2& m){
  output << m.vertices.size() << " ";

  for(unsigned int i = 0; i < m.vertices.size(); i++){
    output << m.vertices[i] << " ";
  }

  output << m.faces.size() << " ";

  for(unsigned int i = 0; i < m.faces.size(); i++){ 
    output << m.faces[i] << " "; 
  } 

  output << m.normals.size() << " ";

  for(unsigned int i = 0; i < m.normals.size(); i++){ 
    output << m.normals[i] << " "; 
  } 

  output << m.face_normals.size() << " ";

  for(unsigned int i = 0; i < m.face_normals.size(); i++){ 
    output << m.face_normals[i] << " "; 
  } 

  return output;
}

inline std::istream &operator>>(std::istream &input, TriangleMesh2& m){

  unsigned int num_vertices;  
  input >> num_vertices;  
  m.vertices.resize(num_vertices);  
  for(unsigned int i = 0; i < num_vertices; i++){
    input >> m.vertices[i];
  }

  unsigned int num_faces;
  input >> num_faces;
  m.faces.resize(num_faces);
  for(unsigned int i = 0; i < num_faces; i++){ 
    input >> m.faces[i]; 
  }   

  unsigned int num_normals;
  input >> num_normals;
  m.normals.resize(num_normals);
  for(unsigned int i = 0; i < num_normals; i++){ 
    input >> m.normals[i]; 
  }   

  unsigned int num_face_normals;
  input >> num_face_normals;
  m.face_normals.resize(num_face_normals);
  for(unsigned int i = 0; i < num_face_normals; i++){ 
    input >> m.face_normals[i]; 
  }   

  return input;
}


void LoadMultipleFromOBJFile(const std::string& filename, std::vector<TriangleMesh2>& meshes, std::vector<std::string>& groupnames);




#endif
