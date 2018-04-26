#ifndef TETMESH_H
#define TETMESH_H

#include "renderable.h"
#include "trianglemesh2.h"
#include <vector>

class IndexQuadruple{
  friend std::ostream &operator<<(std::ostream &output, const IndexQuadruple&);
  friend std::istream &operator>>(std::istream &input,        IndexQuadruple&);
 private:
  Index quadruple[4];
 public:
  IndexQuadruple() { quadruple[0] = 0; quadruple[1] = 0; quadruple[2] = 0; quadruple[3] = 0;}
  IndexQuadruple(Index a,Index b,Index c,Index d) { quadruple[0] = a; quadruple[1] = b; quadruple[2] = c; quadruple[3] = d;}
  
  const Index& A() const {return quadruple[0];}
        Index& A()       {return quadruple[0];}
  const Index& B() const {return quadruple[1];}
        Index& B()       {return quadruple[1];}
  const Index& C() const {return quadruple[2];}
        Index& C()       {return quadruple[2];}
  const Index& D() const {return quadruple[3];}
        Index& D()       {return quadruple[3];}

  const Index& operator[](unsigned int i) const{ return quadruple[i]; }
        Index& operator[](unsigned int i)      { return quadruple[i]; }
  
};

inline std::ostream &operator<<(std::ostream &output, const IndexQuadruple& t){
  output << t.A() << " " << t.B() << " " << t.C() << " " << t.D();
  return output;
}

inline std::istream &operator>>(std::istream &input,        IndexQuadruple& t){
  input >> t.A() >> t.B() >> t.C() >> t.D();
  return input;
}



class SignedIndexQuadruple{
  friend std::ostream &operator<<(std::ostream &output, const SignedIndexQuadruple&);
  friend std::istream &operator>>(std::istream &input,        SignedIndexQuadruple&);
 private:
  SignedIndex quadruple[4];
 public:
  SignedIndexQuadruple() { quadruple[0] = 0; quadruple[1] = 0; quadruple[2] = 0; quadruple[3] = 0;}
  SignedIndexQuadruple(SignedIndex a,SignedIndex b,SignedIndex c,SignedIndex d) { quadruple[0] = a; quadruple[1] = b; quadruple[2] = c; quadruple[3] = d;}
  
  const SignedIndex& A() const {return quadruple[0];}
        SignedIndex& A()       {return quadruple[0];}
  const SignedIndex& B() const {return quadruple[1];}
        SignedIndex& B()       {return quadruple[1];}
  const SignedIndex& C() const {return quadruple[2];}
        SignedIndex& C()       {return quadruple[2];}
  const SignedIndex& D() const {return quadruple[3];}
        SignedIndex& D()       {return quadruple[3];}

  const SignedIndex& operator[](unsigned int i) const{ return quadruple[i]; }
        SignedIndex& operator[](unsigned int i)      { return quadruple[i]; }
  
};

inline std::ostream &operator<<(std::ostream &output, const SignedIndexQuadruple& t){
  output << t.A() << " " << t.B() << " " << t.C() << " " << t.D();
  return output;
}

inline std::istream &operator>>(std::istream &input,        SignedIndexQuadruple& t){
  input >> t.A() >> t.B() >> t.C() >> t.D();
  return input;
}




/*
 * The face opposite vertex 0 is (1,2,3), opposite vertex 1 is (0,2,3), ...
 */
static const unsigned int tet_face_order[4][3] = { {1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}};



class TetMesh{
 public:

  std::vector< Vector3d >             vertices;
  std::vector< IndexTriple >          faces;
  std::vector< IndexQuadruple >       tetrahedra;
  std::vector< SignedIndexQuadruple > neighbors;

  Index num_surfacevertices;

  void LoadFromFile(const std::string& filename);
  void VerifyMesh();
  void ReorderVertices();
  void ComputeNeighbors();
  void LoadNetGenMesh(const std::string& model_name);
  void LoadTetGenMesh(const std::string& model_name);


  
};


#endif
