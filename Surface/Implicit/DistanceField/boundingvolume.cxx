#include "boundingvolume.h"

//#include "openGLutils.h"  

/*
 * Expands this AABB to include a given vector
 */
void AABoundingBox::include(const Vector3d& v){
  if(!initialized){
    v0 = v;
    v1 = v;
    initialized = true;
  } else {
    for(int i = 0; i < 3; i++){
      v0[i] = std::min(v0[i],v[i]);
      v1[i] = std::max(v1[i],v[i]);
    }
  }
}    

/*
 * Expands this BoundingSphere to include a given vector
 */
void BoundingSphere::include(const Vector3d& v){
  if(!initialized){
    center = v;
    radius = 0;
    initialized = true;
  } else {
    scalar r = (v - center).norm();
    if(r > radius){  //point is currently outside
      scalar c = (radius / r);

      center = ((1 + c) * center + (1 - c) * v)/2;
      radius = (radius + r)/2;
    }
  }
}    

/*
 * Expands this AABB to include another AABB
 */
void AABoundingBox::include(const AABoundingBox& b){
  if(!initialized){
    (*this) = b;
    return;
  } else if(!b.initialized){
    return;
  }
  
  for(int i = 0; i < 3; i++){
    v0[i] = std::min(v0[i],b.v0[i]);
    v1[i] = std::max(v1[i],b.v1[i]);
  }
  
}

/*
 * Expands this BoundingSphere to include another BoundingSphere
 */
void BoundingSphere::include(const BoundingSphere& b){
  if(!initialized){
    (*this) = b;
    return;
  } else if(!b.initialized){
    return;
  }
  
  scalar r = (b.center - center).norm();

  if(r + b.radius > radius){
    scalar c1 = (radius / r);
    scalar c2 = (b.radius / r);
    
    center = ((1 + c1 - c2) * center + (1 - c1 + c2) * b.center)/2;
    radius = (radius + r + b.radius)/2;
  }
  
}


void AABoundingBox::Render(){
  glPushMatrix();
  {
    glTranslated(0.5 * (v1.X() - v0.X()) + v0.X(),
		 0.5 * (v1.Y() - v0.Y()) + v0.Y(),
		 0.5 * (v1.Z() - v0.Z()) + v0.Z());
    glScaled(v1.X() - v0.X(),
	     v1.Y() - v0.Y(),
	     v1.Z() - v0.Z());	     
   // glutSolidCube(1.0);
  }
  glPopMatrix();
}

void BoundingSphere::Render(){
  glPushMatrix();
  {
    glTranslated(center.X(), center.Y(), center.Z());
   // glutSolidSphere(radius,8,8);
  }
  glPopMatrix();
}


