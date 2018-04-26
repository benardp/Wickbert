#include "trianglemesh2.h"

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "ctype.h"

#include "defs.h"
#include "utils.h"
#include "opengl.h"
#include "stopwatch.h"

//added because I took it from the .h
#include "matrix3d.h" 
#include "trianglenearestpoint.h"

TriangleMesh2::TriangleMesh2(const std::string filename){
	loadedFlag = LoadFromOBJFile(filename);
	if(!loadedFlag) 
		std::cout << "unable to load mesh" << std::endl;
}


TriangleMesh2::~TriangleMesh2()
{
	loadedFlag = false;
}


Triangle TriangleMesh2::operator[](unsigned int i) const{
  const IndexTriple& t = faces[i];
//   if(faces.size() == face_normals.size()){
//     const IndexTriple& n = face_normals[i];
//     return Triangle(vertices[t.A()],vertices[t.B()],vertices[t.C()],normals[n.A()],normals[n.B()],normals[n.C()]);
//   } else {
    return Triangle(vertices[t.A()],vertices[t.B()],vertices[t.C()]);
    //  }
}


/*
 * Compute the volume of this mesh, assuming the mesh actually defines a volume
 */
scalar TriangleMesh2::Volume() const{
  scalar volume = 0.0;
  
  for(unsigned int i = 0; i < size(); i++){
    const Triangle& t = (*this)[i];
    scalar x1 = t[0].X();
    scalar x2 = t[1].X();
    scalar x3 = t[2].X();
    
    scalar y1 = t[0].Y();
    scalar y2 = t[1].Y();
    scalar y3 = t[2].Y();
    
    scalar z1 = t[0].Z();
    scalar z2 = t[1].Z();
    scalar z3 = t[2].Z();
    volume += ((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)) * (z1 + z2 + z3);
  }					    

  return volume / 6.0;
}


void Subexpressions(scalar w0,scalar w1,scalar w2,scalar& f1,scalar& f2,scalar& f3,scalar& g0,scalar& g1,scalar& g2)
{
  scalar temp0 = w0+w1; f1 = temp0+w2; scalar temp1 = w0*w0; scalar temp2 = temp1+w1*temp0;
  f2 = temp2+w2*f1; f3 = w0*temp1+w1*temp2+w2*f2;
  g0 = f2+w0*(f1+w0); g1 = f2+w1*(f1+w1); g2 = f2+w2*(f1+w2);
}
/* 
 * Determine volume, center of mass, and inertia tensor from a triangle mesh
 *
 * Adapted from code in Dave Eberly's paper "Polyhedral Mass Properties (Revisited)"
 * which is based on Brian Mirtich's paper:
 *    "Fast and Accurate computation of polyhedral mass properties"
 *     Journal of Graphics Tools, vol. 1, no. 2, pp. 31 50, 1996.
 */
void TriangleMesh2::MassProperties (scalar& volume, Vector3d& cm, Matrix3d& inertia) const {
  const scalar mult[10] = {1.0/6.0, 1.0/24.0, 1.0/24.0, 1.0/24.0, 1.0/60.0, 1.0/60.0, 1.0/60.0, 1.0/120.0, 1.0/120.0, 1.0/120.0};
  scalar intg[10] = {0,0,0,0,0,0,0,0,0,0}; // order: 1, x, y, z, x^2, y^2, z^2, xy, yz, zx

  IndexTriple tri;
  Vector3d v0, v1, v2;
  scalar x0, x1, x2, y0, y1, y2, z0, z1, z2;
  scalar a1, b1, c1, a2, b2, c2, d0, d1, d2;
  scalar f1x, f2x, f3x, g0x, g1x, g2x;
  scalar f1y, f2y, f3y, g0y, g1y, g2y;
  scalar f1z, f2z, f3z, g0z, g1z, g2z;

  for (unsigned int t = 0; t < faces.size(); t++){
    // get vertices of triangle t
    tri = faces[t];
    v0 = vertices[tri.A()]; v1 = vertices[tri.B()]; v2 = vertices[tri.C()]; 
    x0 = v0.X(); y0 = v0.Y(); z0 = v0.Z();
    x1 = v1.X(); y1 = v1.Y(); z1 = v1.Z();
    x2 = v2.X(); y2 = v2.Y(); z2 = v2.Z();
    // get edges and cross product of edges
    a1 = x1-x0; b1 = y1-y0; c1 = z1-z0; a2 = x2-x0; b2 = y2-y0; c2 = z2-z0;
    d0 = b1*c2-b2*c1; d1 = a2*c1-a1*c2; d2 = a1*b2-a2*b1;
    // compute integral terms
    Subexpressions(x0,x1,x2,f1x,f2x,f3x,g0x,g1x,g2x);
    Subexpressions(y0,y1,y2,f1y,f2y,f3y,g0y,g1y,g2y);
    Subexpressions(z0,z1,z2,f1z,f2z,f3z,g0z,g1z,g2z);
    // update integrals
    intg[0] += d0*f1x;
    intg[1] += d0*f2x; intg[2] += d1*f2y; intg[3] += d2*f2z;
    intg[4] += d0*f3x; intg[5] += d1*f3y; intg[6] += d2*f3z;
    intg[7] += d0*(y0*g0x+y1*g1x+y2*g2x);
    intg[8] += d1*(z0*g0y+z1*g1y+z2*g2y);
    intg[9] += d2*(x0*g0z+x1*g1z+x2*g2z);
  }
  for (unsigned int i = 0; i < 10; i++)
    intg[i] *= mult[i];

  //  cout << "MULT[0]" << mult[0] << endl << endl;

  // volume
  volume = intg[0];

  
  // center of mass
  cm = Vector3d(intg[1],intg[2],intg[3]) / volume;
  
  // inertia tensor relative to center of mass
  inertia[0][0] = intg[5]+intg[6]-volume*(cm.Y()*cm.Y()+cm.Z()*cm.Z());
  inertia[1][1] = intg[4]+intg[6]-volume*(cm.Z()*cm.Z()+cm.X()*cm.X());
  inertia[2][2] = intg[4]+intg[5]-volume*(cm.X()*cm.X()+cm.Y()*cm.Y());
  inertia[0][1] = -(intg[7]-volume*cm.X()*cm.Y());
  inertia[1][0] = inertia[0][1];
  inertia[1][2] = -(intg[8]-volume*cm.Y()*cm.Z());
  inertia[2][1] = inertia[1][2];
  inertia[0][2] = -(intg[9]-volume*cm.Z()*cm.X());
  inertia[2][0] = inertia[0][2];
}




TriangleMesh2 TriangleMesh2::UnitCube(){
  TriangleMesh2 mesh;
  
  Vector3d nnn(-0.5,-0.5,-0.5);
  Vector3d nnp(-0.5,-0.5, 0.5);
  Vector3d npn(-0.5, 0.5,-0.5);
  Vector3d npp(-0.5, 0.5, 0.5);
  Vector3d pnn( 0.5,-0.5,-0.5);
  Vector3d pnp( 0.5,-0.5, 0.5);
  Vector3d ppn( 0.5, 0.5,-0.5);
  Vector3d ppp( 0.5, 0.5, 0.5);

  mesh.vertices.push_back(Vector3d(-0.5,-0.5, 0.5));
  mesh.vertices.push_back(Vector3d( 0.5,-0.5, 0.5));
  mesh.vertices.push_back(Vector3d(-0.5, 0.5, 0.5));
  mesh.vertices.push_back(Vector3d( 0.5, 0.5, 0.5));
  mesh.vertices.push_back(Vector3d(-0.5, 0.5,-0.5));
  mesh.vertices.push_back(Vector3d( 0.5, 0.5,-0.5));
  mesh.vertices.push_back(Vector3d(-0.5,-0.5,-0.5));
  mesh.vertices.push_back(Vector3d( 0.5,-0.5,-0.5));
  
  mesh.faces.push_back(IndexTriple(0,1,2));
  mesh.faces.push_back(IndexTriple(1,3,2));
  mesh.faces.push_back(IndexTriple(2,3,4));
  mesh.faces.push_back(IndexTriple(3,5,4));
  mesh.faces.push_back(IndexTriple(4,5,6));
  mesh.faces.push_back(IndexTriple(5,7,6));
  mesh.faces.push_back(IndexTriple(6,7,0));
  mesh.faces.push_back(IndexTriple(7,1,0));
  mesh.faces.push_back(IndexTriple(1,7,3));
  mesh.faces.push_back(IndexTriple(7,5,3));
  mesh.faces.push_back(IndexTriple(6,0,4));
  mesh.faces.push_back(IndexTriple(0,2,4));

  mesh.normals.push_back(Vector3d( 0.0, 0.0, 1.0));
  mesh.normals.push_back(Vector3d( 0.0, 0.0,-1.0));
  mesh.normals.push_back(Vector3d( 0.0, 1.0, 0.0));
  mesh.normals.push_back(Vector3d( 0.0,-1.0, 0.0));
  mesh.normals.push_back(Vector3d( 1.0, 0.0, 0.0));
  mesh.normals.push_back(Vector3d(-1.0, 0.0, 0.0));

  mesh.face_normals.push_back(IndexTriple(0,0,0));
  mesh.face_normals.push_back(IndexTriple(0,0,0));
  mesh.face_normals.push_back(IndexTriple(2,2,2));
  mesh.face_normals.push_back(IndexTriple(2,2,2));
  mesh.face_normals.push_back(IndexTriple(1,1,1));
  mesh.face_normals.push_back(IndexTriple(1,1,1));
  mesh.face_normals.push_back(IndexTriple(3,3,3));
  mesh.face_normals.push_back(IndexTriple(3,3,3));
  mesh.face_normals.push_back(IndexTriple(4,4,4));
  mesh.face_normals.push_back(IndexTriple(4,4,4));
  mesh.face_normals.push_back(IndexTriple(5,5,5));
  mesh.face_normals.push_back(IndexTriple(5,5,5));

  return mesh;
}

void TriangleMesh2::Translate(const Vector3d& v){
  for(unsigned int i = 0; i < vertices.size(); i++){
    vertices[i] += v;
  }
}

void TriangleMesh2::Scale(const Vector3d& v){
  for(unsigned int i = 0; i < vertices.size(); i++){
    vertices[i][0] *= v[0];
    vertices[i][1] *= v[1];
    vertices[i][2] *= v[2];
  }
}


AABoundingBox TriangleMesh2::BoundingBox() const{
  AABoundingBox bbox;

  for(unsigned int i = 0; i < vertices.size(); i++){
    bbox.include(vertices[i]);
  }

  return bbox;
}





void LoadMultipleFromOBJFile(const std::string& filename, std::vector<TriangleMesh2>& meshes, std::vector<std::string>& groupnames){
  TriangleMesh2 current;

  std::ifstream file(filename.c_str());
  if(!file){
    std::cerr << "Unable to open OBJ file: " << filename << std::endl;
    exit(EXIT_FAILURE);
  }

  // Remember OBJ indices start at 1 instead of 0!

  bool group_started = false;

  int base_vertex_index = 1;
  int base_normal_index = 1;

  while(!file.eof()){
    std::string str;
    file >> str;

    if(str == "v"){         //vertex 
      if(group_started){
	group_started = false;
	meshes.push_back(current);
	base_vertex_index += current.vertices.size();
	base_normal_index += current.normals.size();
	current.vertices.clear();
	current.normals.clear();
	current.faces.clear();
	current.face_normals.clear();
      }
      Vector3d temp;
      file >> temp;
      current.vertices.push_back(temp);
    }  
    else if(str == "vn"){    //vertex normal
      Vector3d temp;
      file >> temp;
      current.normals.push_back(temp);
    }
    else if(str == "f"){    //face (we only handle triangles)
      bool has_normals = false;
      int v0, v1, v2;
      int vn0, vn1, vn2;

      vn0 = vn1 = vn2 = 0;
      
      file >> v0;
      while(file.peek() == '/'){
	has_normals = true;
	file.ignore(1);
	if(file.peek() != '/')
	  file >> vn0;
      }
      
      file >> v1;
      while((char) file.peek() == '/'){
	file.ignore(1);
	if(file.peek() != '/')
	  file >> vn1;
      }
      
      file >> v2;
      while((char)file.peek() == '/'){
	file.ignore(1);
	if(file.peek() != '/')
	  file >> vn2;
      }


      current.faces.push_back(IndexTriple(v0-base_vertex_index,v1-base_vertex_index,v2-base_vertex_index)); //Fix indices
      if(has_normals){
	current.face_normals.push_back(IndexTriple(vn0 - base_normal_index, vn1 - base_normal_index, vn2 - base_normal_index));
      }
    }
    else if(str == "g"){

      group_started = true;
      std::string groupname;
      file >> groupname;
      std::cout << " GROUP FOUND " << groupname << std::endl;
      groupnames.push_back(groupname);
    }
    else{
      move_to_next_line(file);
    }
  }
  file.close();

  meshes.push_back(current);
}


  
void TriangleMesh2::LoadFromSMESHFile(const std::string& filename){

  std::ifstream smesh_file(filename.c_str());
  if(!smesh_file){ std::cerr << "unable to open mesh file: " << filename << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string line;
  std::istringstream linestream;  


  //Skip comments
  do{
    getline(smesh_file,line);
  } while (line[0] == '#');

  linestream.str(line);
  unsigned int num_vertices;
  linestream >> num_vertices;

#ifdef DEBUG
  std::cout << "num_vertices " << num_vertices << "\n";
#endif

  vertices.resize(num_vertices);
  for(unsigned int i = 0; i < num_vertices; i++){
    getline(smesh_file,line);
    linestream.str(line);
    linestream.clear();

    unsigned int index; linestream >> index; assert(i == index);
        
    linestream >> vertices[i];
  }
  

  //Skip comments
  do{
    line.clear();
    getline(smesh_file,line);
  } while (line[0] == '#');

  
  linestream.str(line);
  linestream.clear();
  unsigned int num_faces;
  linestream >> num_faces;

#ifdef DEBUG
  std::cout << "num faces " << num_faces << "\n";
#endif

  faces.resize(num_faces);
  for(unsigned int i = 0; i < num_faces; i++){
    getline(smesh_file,line);
    linestream.str(line);
    linestream.clear();

    unsigned int num_indices;  linestream >> num_indices;  assert(num_indices == 3);

    linestream >> faces[i];
#ifdef DEBUG
    std::cout << faces[i] << "\n";
#endif
  }

  smesh_file.close();
}


void TriangleMesh2::SaveToOBJFile(const std::string& filename){
  std::ofstream file(filename.c_str());
  if(!file){
    std::cerr << "Unable to open OBJ file: " << filename << std::endl;
    exit(EXIT_FAILURE);
  }

  for(unsigned int i = 0; i < vertices.size(); i++){
    file << "v " << vertices[i] << "\n";
  }

  for(unsigned int i = 0; i < faces.size(); i++){
    file << "f " << (faces[i][0]+1) << " " << (faces[i][1]+1) << " " << (faces[i][2]+1) << "\n";
  }


  file.close();
}

bool TriangleMesh2::LoadFromOBJFile(const std::string& filename){
  std::ifstream file(filename.c_str());
  if(!file){
    std::cerr << "Unable to open OBJ file: " << filename << std::endl;
    
	 return false;
  }

  // Remember OBJ indices start at 1 instead of 0!

  while(!file.eof()){
    std::string str;
    file >> str;

    if(str == "v"){         //vertex 
      Vector3d temp;
      file >> temp;
      vertices.push_back(temp);
    }  
    else if(str == "vn"){    //vertex normal
      Vector3d temp;
      file >> temp;
      normals.push_back(temp);
    }
    else if(str == "f"){    //face (we only handle triangles)
      bool has_normals = false;
      int v0, v1, v2;
      int vn0, vn1, vn2;

      vn0 = vn1 = vn2 = 0;
      
      file >> v0;
      while(file.peek() == '/'){
	has_normals = true;
	file.ignore(1);
	if(file.peek() != '/')
	  file >> vn0;
      }
      
      file >> v1;
      while((char) file.peek() == '/'){
	file.ignore(1);
	if(file.peek() != '/')
	  file >> vn1;
      }
      
      file >> v2;
      while((char)file.peek() == '/'){
	file.ignore(1);
	if(file.peek() != '/')
	  file >> vn2;
      }

      faces.push_back(IndexTriple(v0-1,v1-1,v2-1)); //Fix indices
      if(has_normals){
	face_normals.push_back(IndexTriple(vn0 - 1, vn1 - 1, vn2 - 1));
      }
    }
    else{
      move_to_next_line(file);
    }
  }
  file.close();

  std::cout << "Loaded " << faces.size() << " triangles from " << filename << std::endl;
  return true;
}


void TriangleMesh2::Render(){
  const TriangleMesh2& mesh = (*this);

  for(unsigned int i = 0; i < mesh.size(); i++){
    mesh[i].Render();
  }
}
