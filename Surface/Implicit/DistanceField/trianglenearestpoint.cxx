#include "trianglenearestpoint.h"


TriangleNearestPoint::TriangleNearestPoint(const Triangle& T) {
	initFromVertices(T.A(),T.B(),T.C());
}	

TriangleNearestPoint::TriangleNearestPoint(const gmVector3& a, const gmVector3& b, const gmVector3& c) {
	initFromVertices	(	Vector3d(a[0],a[1],a[2]),
							Vector3d(b[0],b[1],b[2]),
							Vector3d(c[1],c[2],c[3])
						);
}	
	
	
void TriangleNearestPoint::initFromVertices(const Vector3d &a, const Vector3d &b, const Vector3d &c)
{
  corners[0] = a;
  corners[1] = b;
  corners[2] = c;

  span01 = (corners[1] - corners[0]);
  span02 = (corners[2] - corners[0]);
  span12 = (corners[2] - corners[1]);

  unitnormal = cross(span01,span02);
  
  if (unitnormal.norm2() > 0.0)
    normalize(unitnormal);
  else
    throw DegenerateTriangleException();
  
  u11 = span01 * span01;
  u22 = span02 * span02;
  u33 = span12 * span12;  

  u12 = span01 * span02;
  u13 = span01 * span12;

  u_span01 = span01 / u11;
  u_span02 = span02 / u22;
  u_span12 = span12 / u33;
  

  scalar d = u11 * u22 - u12 * u12;
  d_inv = (scalar) 1.0/d;
}

Vector3d TriangleNearestPoint::NearestPointTo(const Vector3d& v) const{
  VoronoiRegion temp;
  return NearestPointToVerbose(v,temp);
}

Vector3d TriangleNearestPoint::NearestPointToVerbose(const Vector3d& v, VoronoiRegion & voronoi_region) const{
  Vector3d p = v - corners[0];

  //project to triangle's plane
  p -= (p*unitnormal) * unitnormal;  

  const scalar p01 = p * span01;
  const scalar p02 = p * span02;
  //  scalar p12 = (v - corners[1] - u*truenormal) * span12;
  const scalar p12 = p*span12 - u13;


  const scalar z1 = d_inv * (p01*u22 - p02*u12);
  const scalar z2 = d_inv * (p02*u11 - p01*u12);
  const scalar z0 = 1 - z1 - z2;
  
  int barycentric_region = 0;
  if(z0 > 0.0) barycentric_region += 1;
  if(z1 > 0.0) barycentric_region += 2;
  if(z2 > 0.0) barycentric_region += 4;

  switch(barycentric_region){
  case 1:
    if(p01 >= u11){
      voronoi_region = CORNER_B_REGION;
      return corners[1];
    } else if (p02 >= u22){
      voronoi_region = CORNER_C_REGION;
      return corners[2];
    } else if (p01 < 0) {
      if(p02 < 0){ 
	voronoi_region = CORNER_A_REGION;
	return corners[0];
      } else {
	voronoi_region = EDGE_AC_REGION;
	return corners[0] + p02 * u_span02;
      }
    } else {
      voronoi_region = EDGE_AB_REGION;
      return corners[0] + p01 * u_span01;
    }    

  case 2:
    if(p01 <= 0){
      voronoi_region = CORNER_A_REGION;
      return corners[0];
    } else if (p12 >= u33){
      voronoi_region = CORNER_C_REGION;
      return corners[2];
    } else if (p01 < u11) {
      voronoi_region = EDGE_AB_REGION;
      return corners[0] + p01 * u_span01;
    } else if (p12 > 0 ) {
      voronoi_region = EDGE_BC_REGION;
      return corners[1] + p12 * u_span12;
    } else {
      voronoi_region = CORNER_B_REGION;
      return corners[1];
    }

  case 3:
    if (p01 <= 0.0){
      voronoi_region = CORNER_A_REGION;
      return corners[0];
    } else if (p01 >= u11){
      voronoi_region = CORNER_B_REGION;
      return corners[1];
    } else {
      voronoi_region = EDGE_AB_REGION;
      return corners[0] + p01 * u_span01;
    }
      
  case 4:
    if(p02 <= 0){
      voronoi_region = CORNER_A_REGION;
      return corners[0];
    } else if (p12 <= 0) {
      voronoi_region = CORNER_B_REGION;
      return corners[1];
    } else if (p02 < u22) {
      voronoi_region = EDGE_AC_REGION;
      return corners[0] + p02 * u_span02;
    } else if (p12 < u33){
      voronoi_region = EDGE_BC_REGION;
      return corners[1] + p12 * u_span12;
    } else {
      voronoi_region = CORNER_C_REGION;
      return corners[2];
    }

  case 5:
    if (p02 <= 0.0){
      voronoi_region = CORNER_A_REGION;
      return corners[0];
    } else if (p02 >= u22){
      voronoi_region = CORNER_C_REGION;
      return corners[2];
    } else {
      voronoi_region = EDGE_AC_REGION;
      return corners[0] + p02 * u_span02;
    }

  case 6:
    if(p12 <= 0){
      voronoi_region = CORNER_B_REGION;
      return corners[1];
    } else{
      if(p12 >= u33){
	voronoi_region = CORNER_C_REGION;
	return corners[2];
      } else{
	voronoi_region = EDGE_BC_REGION;
	return corners[1] + p12 * u_span12;
      }
    }
  
  case 7: //inside triangle
    voronoi_region = FACE_REGION;
    return corners[0] + p;  // same as: corners[0] + z1*span01 + z2*span02;

  default:
    assert(false);
    return Vector3d(0,0,0);
  }
}


gmVector3 TriangleNearestPoint::NearestPointTo(const gmVector3 & pos, VoronoiRegion & region) const
{
	Vector3d p(pos[0],pos[1],pos[2]);
	Vector3d result=NearestPointToVerbose(p,region);
	return gmVector3(result[0],result[1],result[2]);
}



bool TriangleNearestPoint::Test(){

  unsigned int NUM_TRIANGLES = 1000;

  const scalar TEST_EPSILON = 1e-5;

  /*
   * Test Nearest Point Method
   */
  std::cout << "TriangleNearestPoint::NearestPointTo(): ";
  for(unsigned int i = 0; i < NUM_TRIANGLES; i++){
    Vector3d A(rand_scalar(2.0) - 1.0,rand_scalar(2.0) - 1.0,rand_scalar(2.0) - 1.0);
    Vector3d B(rand_scalar(2.0) - 1.0,rand_scalar(2.0) - 1.0,rand_scalar(2.0) - 1.0);
    Vector3d C(rand_scalar(2.0) - 1.0,rand_scalar(2.0) - 1.0,rand_scalar(2.0) - 1.0);

    Triangle T(A,B,C);
    TriangleNearestPoint TNP(T);

    unsigned int NUM_SAMPLES = 1000;
    
    for(unsigned int j = 0; j < NUM_SAMPLES; j++){
      Vector3d sample(rand_scalar(3.0) - 1.5,rand_scalar(3.0) - 1.5,rand_scalar(3.0) - 1.5);

      Vector3d nearest_pt = TNP.NearestPointTo(sample);
      scalar distance_squared = (sample-nearest_pt).norm2();

      //Make sure corners aren't any closer than nearest pt   
      assert((sample - T.A()).norm2() > distance_squared - TEST_EPSILON);
      assert((sample - T.B()).norm2() > distance_squared - TEST_EPSILON);
      assert((sample - T.C()).norm2() > distance_squared - TEST_EPSILON);

      //Make sure edges aren't any closer than nearest pt   
      Vector3d edge01 = T.B() - T.A();
      Vector3d edge12 = T.C() - T.B();
      Vector3d edge20 = T.A() - T.C();
      
      Vector3d s0 = sample - T.A();
      Vector3d s1 = sample - T.B();
      Vector3d s2 = sample - T.C();

      scalar p0 = edge01 * s0;
      scalar p1 = edge12 * s1;
      scalar p2 = edge20 * s2;

      if(p0 >= 0 && p0 <= edge01.norm2()){
	Vector3d n0 = T.A() + (p0 / edge01.norm2()) * edge01;
	assert((sample - n0).norm2() > distance_squared - TEST_EPSILON);
      }
      if(p1 >= 0 && p1 <= edge12.norm2()){
	Vector3d n1 = T.B() + (p1 / edge12.norm2()) * edge12;
	assert((sample - n1).norm2() > distance_squared - TEST_EPSILON);
      }
      if(p2 >= 0 && p2 <= edge20.norm2()){
	Vector3d n2 = T.C() + (p2 / edge20.norm2()) * edge20;
	assert((sample - n2).norm2() > distance_squared - TEST_EPSILON);
      }
    }

    //Now test interior points
    for(unsigned int j = 0; j < NUM_SAMPLES; j++){
      Vector3d point_in_triangle = T.RandomSample();
      Vector3d sample = point_in_triangle + (rand_scalar(2.0) - 1.0) * T.Normal();
      
      Vector3d nearest = TNP.NearestPointTo(sample);
      assert((nearest - point_in_triangle).norm2() < TEST_EPSILON);
    }
  }

  std::cout << " passed" << std::endl;

  return true;
}
