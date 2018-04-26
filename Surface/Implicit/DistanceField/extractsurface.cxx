#include "extractsurface.h"


/* 
 *  Shared Methods
 *
 */
void TriangulateCell(const Vector3d * corners, const scalar * values, 
		     const scalar isosurface, 
		     std::vector< Triangle >& triangles);

void TriangulateTetrahedron(const Vector3d& v0, const Vector3d& v1,
			    const Vector3d& v2, const Vector3d& v3,
			    const scalar& val0, const scalar& val1,
			    const scalar& val2, const scalar& val3,
			    const scalar isosurface, 
			    std::vector< Triangle >& triangles);

Vector3d InterpolateVertex(const scalar iso, const Vector3d& p1, const Vector3d& p2, 
			   const scalar valp1, const scalar valp2);

/*
 * ADF Specific Methods
 *
 */
void TriangulateOctree(const OctreeNode * node, scalar isosurface, std::vector< Triangle >& triangles);
void TriangulateOctreeCell(const OctreeNode * node, scalar isosurface, std::vector< Triangle >& triangles);
void ExtractIsosurface(const ADF& adf, scalar isosurface, std::vector< Triangle >& triangles){
  for(ADF::const_iterator i = adf.begin(); i != adf.end(); i++){
    TriangulateOctree(i->second,isosurface,triangles);
  }
}

void TriangulateOctree(const OctreeNode * node, scalar isosurface, std::vector< Triangle >& triangles){
  if(node->is_leaf){
    TriangulateOctreeCell(node,isosurface,triangles);
  } else {
    for(int i = 0; i < 8; i++){
      TriangulateOctree(node->children[i],isosurface,triangles);
    }
  }
}


/*
 * Subdivide this cube into 6 tetrahedra to be triangulated
 */
void TriangulateOctreeCell(const OctreeNode * node, scalar isosurface, std::vector< Triangle >& triangles){
  const Vector3d& corner = node->corner;
  const scalar& edge_length = node->edge_length;

  //Initialize corners of cell
  Vector3d corners[8];
  for(int i = 0; i < 8; i++){ corners[i] = corner; }
  corners[1][2] += edge_length;  corners[2][1] += edge_length;  corners[4][0] += edge_length;
  corners[3][2] += edge_length;  corners[3][1] += edge_length;  corners[5][0] += edge_length;
  corners[5][2] += edge_length;  corners[6][1] += edge_length;  corners[6][0] += edge_length;
  corners[7][2] += edge_length;  corners[7][1] += edge_length;  corners[7][0] += edge_length;
  
  TriangulateCell(corners,node->distances,isosurface,triangles); 
}





/*
 * SparseScalarLattice Specific Methods
 *
 */
void TriangulateLatticeCell(const SparseScalarLattice& lattice, const Cell& cell, scalar isosurface,
			    std::vector< Triangle >& triangles);
void ExtractIsosurface(const SparseScalarLattice& lattice, scalar isosurface, std::vector< Triangle >& triangles){
  CellSet completed_cells;

  SparseScalarLattice::const_iterator i;
  SparseScalarLattice::const_iterator end = lattice.end();
  for(i = lattice.begin(); i != end; i++){
    TriangulateLatticeCell(lattice, i->first, isosurface, triangles);
  }
}

void TriangulateLatticeCell(const SparseScalarLattice& lattice, const Cell& cell, scalar isosurface,
			    std::vector< Triangle >& triangles){

  //Initialize distances, bail out if any corner is missing
  scalar distances[8];
  int n = 0;
  for(int i = 0; i <= 1; i++){
    for(int j = 0; j <= 1; j++){
      for(int k = 0; k <= 1; k++){
	Cell current = cell + Cell(i,j,k);
	if(!lattice.HasValue(current)){
	  return;
	} else {
	  distances[n++] = lattice.Value(current);
	}
      }
    }
  }

  const scalar edge_length = lattice.VoxelSize();
  const Vector3d corner = edge_length * cell.toVector3d();

  //Initialize corners of cell
  Vector3d corners[8];
  for(int i = 0; i < 8; i++){ corners[i] = corner; }
  corners[1][2] += edge_length;  corners[2][1] += edge_length;  corners[4][0] += edge_length;
  corners[3][2] += edge_length;  corners[3][1] += edge_length;  corners[5][0] += edge_length;
  corners[5][2] += edge_length;  corners[6][1] += edge_length;  corners[6][0] += edge_length;
  corners[7][2] += edge_length;  corners[7][1] += edge_length;  corners[7][0] += edge_length;

  TriangulateCell(corners,distances,isosurface,triangles); 
}











/*
 * ---  Shared Methods  ---
 *
 */




/* 
 *  Split a cube into 6 tetrahedra
 */
void TriangulateCell(const Vector3d * corners, const scalar * values, 
		     const scalar isosurface, 
		     std::vector< Triangle >& triangles){

  //Green
  TriangulateTetrahedron(corners[0],corners[1],corners[5],corners[2],
			 values[0],values[1],values[5],values[2],
			 isosurface,triangles);
  TriangulateTetrahedron(corners[0],corners[4],corners[2],corners[5],
			 values[0],values[4],values[2],values[5],
			 isosurface,triangles);
  //Blue
  TriangulateTetrahedron(corners[2],corners[1],corners[5],corners[3],
			 values[2],values[1],values[5],values[3],
			 isosurface,triangles);
  TriangulateTetrahedron(corners[2],corners[4],corners[6],corners[5],
			 values[2],values[4],values[6],values[5],
			 isosurface,triangles);
  //Red
  TriangulateTetrahedron(corners[2],corners[7],corners[3],corners[5],
			 values[2],values[7],values[3],values[5],
			 isosurface,triangles);
  TriangulateTetrahedron(corners[2],corners[6],corners[7],corners[5],
			 values[2],values[6],values[7],values[5],
			 isosurface,triangles);
}

/*
 * Use linear interpolation to approximate the point on the line between p1 and p2
 * which has a given value (iso)
 */
Vector3d InterpolateVertices(const scalar iso, const Vector3d& p1, const Vector3d& p2, scalar valp1, scalar valp2){ 
  scalar mu = (iso - valp1) / (valp2 - valp1);

  return p1 + mu * (p2 - p1);
}


/*
 *  Apply Marching Tetrahedra to a single tetrahedron
 *
 *  Triangles are oriented counter-clockwise
 */
void TriangulateTetrahedron(const Vector3d& v0, const Vector3d& v1,
			    const Vector3d& v2, const Vector3d& v3,
			    const scalar& val0, const scalar& val1,
			    const scalar& val2, const scalar& val3,
			    const scalar iso, std::vector< Triangle >& triangles){
  
  int triindex = 0;
  int num_triangles = 0;
 

  Vector3d triangle_vertices[2][3];

  if (val0 < iso) triindex |= 1;
  if (val1 < iso) triindex |= 2;
  if (val2 < iso) triindex |= 4;
  if (val3 < iso) triindex |= 8;

  switch (triindex) {
    //Type 1
  case 0x00:
  case 0x0F:
    num_triangles = 0;
    break;
    
    //Type 2
  case 0x0E:
    num_triangles = 1;
    triangle_vertices[0][0] = InterpolateVertices(iso,v0,v2,val0,val2);
    triangle_vertices[0][1] = InterpolateVertices(iso,v0,v1,val0,val1);				  				  
    triangle_vertices[0][2] = InterpolateVertices(iso,v0,v3,val0,val3);
    break;
    
  case 0x01:
    num_triangles = 1;
    triangle_vertices[0][0] = InterpolateVertices(iso,v0,v1,val0,val1);				  
    triangle_vertices[0][1] = InterpolateVertices(iso,v0,v2,val0,val2);
    triangle_vertices[0][2] = InterpolateVertices(iso,v0,v3,val0,val3);
    break;
    
    
    //Type 3
  case 0x0D:
    num_triangles = 1;
    triangle_vertices[0][0] = InterpolateVertices(iso,v1,v3,val1,val3);
    triangle_vertices[0][1] = InterpolateVertices(iso,v1,v0,val1,val0);
    triangle_vertices[0][2] = InterpolateVertices(iso,v1,v2,val1,val2);				 
    break;
    
  case 0x02:
    num_triangles = 1;
    triangle_vertices[0][0] = InterpolateVertices(iso,v1,v0,val1,val0);
    triangle_vertices[0][1] = InterpolateVertices(iso,v1,v3,val1,val3);
    triangle_vertices[0][2] = InterpolateVertices(iso,v1,v2,val1,val2);
    break;
    
    
    //Type 4
  case 0x0C:
    num_triangles = 2;
    triangle_vertices[0][0] = InterpolateVertices(iso,v0,v3,val0,val3);
    triangle_vertices[0][1] = InterpolateVertices(iso,v0,v2,val0,val2);
    triangle_vertices[0][2] = InterpolateVertices(iso,v1,v3,val1,val3);
    
    triangle_vertices[1][0] = InterpolateVertices(iso,v1,v2,val1,val2);
    triangle_vertices[1][1] = InterpolateVertices(iso,v1,v3,val1,val3);
    triangle_vertices[1][2] = InterpolateVertices(iso,v0,v2,val0,val2);
    break;
    
  case 0x03:
    num_triangles = 2;
    triangle_vertices[0][0] = InterpolateVertices(iso,v0,v2,val0,val2);
    triangle_vertices[0][1] = InterpolateVertices(iso,v0,v3,val0,val3);
    triangle_vertices[0][2] = InterpolateVertices(iso,v1,v3,val1,val3);
    
    triangle_vertices[1][0] = InterpolateVertices(iso,v1,v3,val1,val3);
    triangle_vertices[1][1] = InterpolateVertices(iso,v1,v2,val1,val2);
    triangle_vertices[1][2] = InterpolateVertices(iso,v0,v2,val0,val2);				   
    break;
    
    
    //Type 5
  case 0x0B:
    num_triangles = 1;
    triangle_vertices[0][0] = InterpolateVertices(iso,v2,v1,val2,val1);
    triangle_vertices[0][1] = InterpolateVertices(iso,v2,v0,val2,val0);
    triangle_vertices[0][2] = InterpolateVertices(iso,v2,v3,val2,val3);
    break;
    
  case 0x04:
    num_triangles = 1;
    triangle_vertices[0][0] = InterpolateVertices(iso,v2,v0,val2,val0);
    triangle_vertices[0][1] = InterpolateVertices(iso,v2,v1,val2,val1);
    triangle_vertices[0][2] = InterpolateVertices(iso,v2,v3,val2,val3); 
    break;
    
    
    //Type 6
  case 0x0A:
    num_triangles = 2;
    triangle_vertices[0][0] = InterpolateVertices(iso,v2,v3,val2,val3);
    triangle_vertices[0][1] = InterpolateVertices(iso,v0,v1,val0,val1);
    triangle_vertices[0][2] = InterpolateVertices(iso,v0,v3,val0,val3);			       
    
    triangle_vertices[1][0] = InterpolateVertices(iso,v1,v2,val1,val2);
    triangle_vertices[1][1] = InterpolateVertices(iso,v0,v1,val0,val1);
    triangle_vertices[1][2] = InterpolateVertices(iso,v2,v3,val2,val3);
    break;
    
  case 0x05:
    num_triangles = 2;
    triangle_vertices[0][0] = InterpolateVertices(iso,v0,v1,val0,val1);
    triangle_vertices[0][1] = InterpolateVertices(iso,v2,v3,val2,val3);
    triangle_vertices[0][2] = InterpolateVertices(iso,v0,v3,val0,val3);
    
    triangle_vertices[1][0] = InterpolateVertices(iso,v0,v1,val0,val1);
    triangle_vertices[1][1] = InterpolateVertices(iso,v1,v2,val1,val2);
    triangle_vertices[1][2] = InterpolateVertices(iso,v2,v3,val2,val3);
    break;
    
    
    //Type 7
  case 0x09:
    num_triangles = 2;
    triangle_vertices[0][0] = InterpolateVertices(iso,v1,v3,val1,val3);
    triangle_vertices[0][1] = InterpolateVertices(iso,v0,v1,val0,val1);
    triangle_vertices[0][2] = InterpolateVertices(iso,v2,v3,val2,val3);
    
    triangle_vertices[1][0] = InterpolateVertices(iso,v0,v1,val0,val1);
    triangle_vertices[1][1] = InterpolateVertices(iso,v0,v2,val0,val2);
    triangle_vertices[1][2] = InterpolateVertices(iso,v2,v3,val2,val3);
    break;
    
  case 0x06:
    num_triangles = 2;
    triangle_vertices[0][0] = InterpolateVertices(iso,v0,v1,val0,val1);
    triangle_vertices[0][1] = InterpolateVertices(iso,v1,v3,val1,val3);
    triangle_vertices[0][2] = InterpolateVertices(iso,v2,v3,val2,val3);
    
    triangle_vertices[1][0] = InterpolateVertices(iso,v0,v2,val0,val2);
    triangle_vertices[1][1] = InterpolateVertices(iso,v0,v1,val0,val1);
    triangle_vertices[1][2] = InterpolateVertices(iso,v2,v3,val2,val3);
    break;
    
    
    //Type 8
  case 0x07:
    num_triangles = 1;
    triangle_vertices[0][0] = InterpolateVertices(iso,v3,v2,val3,val2);
    triangle_vertices[0][1] = InterpolateVertices(iso,v3,v0,val3,val0);
    triangle_vertices[0][2] = InterpolateVertices(iso,v3,v1,val3,val1);
    break;
    
  case 0x08:
    num_triangles = 1;
    triangle_vertices[0][0] = InterpolateVertices(iso,v3,v0,val3,val0);
    triangle_vertices[0][1] = InterpolateVertices(iso,v3,v2,val3,val2);
    triangle_vertices[0][2] = InterpolateVertices(iso,v3,v1,val3,val1);
    break;
  }
  

  for(int i = 0; i < num_triangles; i++){
    triangles.push_back(Triangle(triangle_vertices[i][0],
				 triangle_vertices[i][1],
				 triangle_vertices[i][2]));
  }
}

