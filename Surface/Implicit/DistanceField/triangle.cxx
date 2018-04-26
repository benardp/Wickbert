#include "triangle.h"

#include "opengl_utils.h"

Triangle::Triangle(const Vector3d& v0, const Vector3d& v1, const Vector3d& v2){
  corners[0] = v0;
  corners[1] = v1;
  corners[2] = v2;
}

//Triangle::Triangle(const Vector3d& v0, const Vector3d& v1, const Vector3d& v2, const Vector3d& n0, const Vector3d& n1,const Vector3d& n2){
//  corners[0] = v0;
//  corners[1] = v1;
//  corners[2] = v2;
//}

Vector3d Triangle::Normal() const{
  return cross(corners[1] - corners[0], corners[2] - corners[0]);
}


/*
 * Return the perimeter of this triangle
 */
scalar Triangle::Perimeter() const{
  scalar edge01len = (corners[1] - corners[0]).norm();
  scalar edge12len = (corners[2] - corners[1]).norm();
  scalar edge20len = (corners[0] - corners[2]).norm();

  return edge01len + edge12len + edge20len;
}

/*
 * Return the surface area of this triangle
 */
scalar Triangle::SurfaceArea() const{
  return 0.5*Normal().norm();
}

/*
 * Pick a random point within this triangle
 */
Vector3d Triangle::RandomSample() const{
  scalar a = rand_scalar(1.0);
  scalar b = rand_scalar(1.0);

  if(a+b > 1.0){
    a = 1 - a;
    b = 1 - b;
  }

  scalar c = 1 - a - b;

  return Bary(a,b,c);
}

/*
 * Return the point in this triangle corresponding to 
 * the Barycentric coordinates (a,b,c)
 *
 */
Vector3d Triangle::Bary(scalar a, scalar b, scalar c) const{
  return a*corners[0] + b*corners[1] + c*corners[2];
}

/* Storage-free triangle intersection/ray method, adapted from code by:
 *   
 *  Tomas Möller and Ben Trumbore. 
 *  Fast, minimum storage ray-triangle intersection. 
 *  Journal of graphics tools, 2(1):21-28, 1997
 *  http://www.acm.org/jgt/papers/MollerTrumbore97/
 *
 */
bool Triangle::Intersection(const Ray& ray, scalar& t) const{
#define EPSILON 0.000001

  Vector3d tvec, pvec, qvec;
  scalar det,inv_det;

  scalar u, v;
      
  /* find vectors for two edges sharing corner 0 */
  Vector3d edge1 = corners[1] - corners[0];
  Vector3d edge2 = corners[2] - corners[0];
  
  /* begin calculating determinant - also used to calculate U parameter */
  pvec =  cross(ray.d, edge2);
  
  /* if determinant is near zero, ray lies in plane of triangle */
  det = edge1 * pvec;
  
  if (det > -EPSILON && det < EPSILON)
    return false;
  inv_det = 1.0 / det;
  
  /* calculate distance from vert0 to ray origin */
  tvec = ray.o - corners[0];
  
  /* calculate U parameter and test bounds */
  u = inv_det * (tvec * pvec);
  if (u < 0.0 || u > 1.0)
    return false;
  
  /* prepare to test V parameter */
  qvec = cross(tvec, edge1);
  
  /* calculate V parameter and test bounds */
  v = inv_det * (ray.d * qvec);
  if (v < 0.0 || u + v > 1.0)
    return false;
  
  /* calculate t, ray intersects triangle */
  t = inv_det * (edge2 * qvec); // <- can be negative
  
  return true;
}


AABoundingBox Triangle::BoundingBox() const{
  AABoundingBox aabb;

  aabb.include(corners[0]);
  aabb.include(corners[1]);
  aabb.include(corners[2]);
  
  return aabb;
}

bool Triangle::Intersects(const Ray& ray) const{
  scalar temp;
  return Intersection(ray,temp) && temp >= 0.0;
}



/*
 * Translate this triangle by given vector
 */
void Triangle::Translate(const Vector3d& v){
  corners[0] += v;
  corners[1] += v;
  corners[2] += v;
}

void Triangle::Render(){
  Vector3d normal = Normal();
  
  if(normal.norm2() > 0)
    normalize(normal);

  glBegin(GL_TRIANGLES);
  {
    glNormal(normal);
    glVertex(corners[0]);
    glNormal(normal);
    glVertex(corners[1]);
    glNormal(normal);
    glVertex(corners[2]);
  }
  glEnd();

//   glPushMatrix(); 
//   {
//     scalar radius = Perimeter() / 10.0;
//     Vector3d loc = Bary(1.0/3.0,1.0/3.0,1.0/3.0) + radius*(Normal().normalized());
//     glTranslate(loc);
//     glutSolidSphere(radius,4,4);
//   }
//   glPopMatrix();
}



