#include "TriMesh.h"

using namespace std;
using namespace OpenMesh;
typedef MyMesh::Point Point;

// default constructor
TriMesh::TriMesh()
{
	// init something
}

bool TriMesh::readFile(const char* filename)
{
	if (filename==NULL)
		return false;
	mesh.clear();
	OpenMesh::IO::read_mesh(mesh, filename);
	mesh.update_face_normals();
	mesh.update_vertex_normals();
	computeFaceCenter();
	computeBoundingBox();
	objectName = filename;
	return true;
}

bool TriMesh::writeFile(const char* filename)
{
	if (filename==NULL)
		return false;
	mesh.clear();
	OpenMesh::IO::write_mesh(mesh, filename);
	return true;	
}

void TriMesh::computeBoundingBox()
{
	MyMesh::VertexIter vit = mesh.vertices_begin(), vend = mesh.vertices_end();

	MyMesh::Point firstPoint = mesh.point(vit);
	bbox_min = bbox_max = firstPoint;

	for ( ; vit != vend; ++vit)
	{		
		MyMesh::Point p = mesh.point(vit);
		for (int i=0; i < 3; i++)
		{
			bbox_min[i] = ((bbox_min[i]<p[i])?bbox_min[i]:p[i]);
			bbox_max[i] = ((bbox_max[i]>p[i])?bbox_max[i]:p[i]);
		}
	}
}

void TriMesh::computeFaceCenter()
{
	for (MyMesh::FaceIter fit = mesh.faces_begin(); fit != mesh.faces_end(); ++fit)
	{
		Point center=Point(0,0,0);
		for (MyMesh::FaceVertexIter fv_it=mesh.fv_iter(fit.handle()); fv_it ; ++fv_it)
		{
			center+=mesh.point(fv_it);
		}
		mesh.data(*fit).center=center/3;
	}
}

Point TriMesh::projectPointToPlane(const Point &p, const Point &a, const Point &b, const Point &c)
{
	Point n=normalOf(a,b,c);
	n.normalize();
	float t=dot(a,n)-dot(p,n);
	return t*n+p;
}

Point TriMesh::projectPointToPlane(const Point &p, const Point &a, const Point &n)
{
	float t=dot(a,n)-dot(p,n);
	return t*n+p;
}

// todo: use face normals
Point TriMesh::normalOf(const Point &a, const Point &b, const Point &c)
{
	Point e1(b-a);
	Point e2(c-a);
	// cross product
	return e1%e2;
}

// test if a point is in triangle
bool TriMesh::isInTriangle(const Point &p, MyMesh::FaceHandle f)
{
	// check for in plane
	MyMesh::FaceVertexIter fv_it=mesh.fv_iter(f);
	MyMesh::VertexHandle va=(fv_it.handle());
	MyMesh::VertexHandle vb=((++fv_it).handle());
	MyMesh::VertexHandle vc=((++fv_it).handle());
	
	const Point &pa=mesh.point(va);
	const Point &pb=mesh.point(vb);
	const Point &pc=mesh.point(vc);
	// make sure p q is in a,b,c
	return isInTriangle(p,pa,pb,pc);
}

bool TriMesh::isInTriangle(MyMesh::VertexHandle p, MyMesh::VertexHandle a, MyMesh::VertexHandle b, MyMesh::VertexHandle c)
{
	const Point &pp=mesh.point(p);
	const Point &pa=mesh.point(a);
	const Point &pb=mesh.point(b);
	const Point &pc=mesh.point(c);
	return isInTriangle(pp,pa,pb,pc);
}

// http://www.blackpawn.com/texts/pointinpoly/default.html
bool TriMesh::isInTriangle(const MyMesh::Point &p, const MyMesh::Point &a, const MyMesh::Point &b, const MyMesh::Point &c)
{
	return (
        isOnSameSide(p,a,b,c) && 
		isOnSameSide(p,b,a,c) && 
		isOnSameSide(p,c,a,b) );
}

bool TriMesh::isOnSameSide(const MyMesh::Point &p1, const MyMesh::Point &p2, const MyMesh::Point &v1, const MyMesh::Point &v2)
{
	Point a=(v2-v1)%(p1-v1);
	Point b=(v2-v1)%(p2-v1);
	return (OpenMesh::dot(a,b) >= 0);
}

// closest distance to a point
double TriMesh::distancePointToMesh(const Point &p)
{
	//MyMesh::FaceHandle fh=closestFace(p);
	//return (p-mesh.face(fh).center).length();

	//MyMesh::VertexHandle vh=closestVertex(p);
	//return (p-mesh.point(vh)).length();

	//double dist=1e50;
	//for (MyMesh::FaceIter fit = mesh.faces_begin(); fit != mesh.faces_end(); ++fit)
	//{
	//	double temp=distancePointToTriangle(p,fit.handle());
	//	if (temp<dist)
	//		dist=temp;
	//}
	//return dist;

	MyMesh::FaceHandle fh=closestFace(p);
	double dist=distancePointToTriangle(p,fh);
	// offset
	return dist;
}

double TriMesh::distancePointToTriangle(const Point &p, MyMesh::FaceHandle fh)
{
	double dist=1e50;

	MyMesh::FaceVertexIter fv_it=mesh.fv_iter(fh);
	MyMesh::VertexHandle va=(fv_it.handle());
	MyMesh::VertexHandle vb=((++fv_it).handle());
	MyMesh::VertexHandle vc=((++fv_it).handle());
	
	const Point &pa=mesh.point(va);
	const Point &pb=mesh.point(vb);
	const Point &pc=mesh.point(vc);

	if (isInTriangle(p,pa,pb,pc))
	{
		Point pq=projectPointToPlane(p,pa,pb,pc);
		dist=(pq-p).length();
	}
	else
	{
		double dist0=distancePointToLine(p,pa,pb);
		double dist1=distancePointToLine(p,pb,pc);
		double dist2=distancePointToLine(p,pc,pa);
		if (dist0<dist) dist=dist0;
		if (dist1<dist) dist=dist1;
		if (dist2<dist) dist=dist2;
	}
	return dist;
}

//gmVector3 TriMesh::grad(const gmVector3 & x)
//{
//	return normal(x);
//}

gmVector3 TriMesh::normal(const gmVector3 & x)
{
	// find closest point;

	MyMesh::FaceHandle fh=closestFace(Point((float)x[0],(float)x[1],(float)x[2]));
	
	Point n=mesh.normal(fh);
	return gmVector3(n[0],n[1],n[2]);
}

// http://www.magic-software.com/Documentation/DistancePointLine.pdf 
float TriMesh::distancePointToLine(const Point &p, const Point &a, const Point &b)
{
	Point ab=b-a;
	float t=OpenMesh::dot((p-a),ab)/ab.sqrnorm();
	if (t<=0)
		return (p-a).length();
	else if (t>=1)
		return (p-b).length();
	else
		return (p-(a+t*ab)).length();
}

MyMesh::FaceHandle TriMesh::closestFace(const Point &p)
{
	// find closest point;
	double distance=1.0e50;
	MyMesh::FaceHandle result;
	for(MyMesh::FaceIter fit = mesh.faces_begin(); fit != mesh.faces_end(); ++fit)
	{
		MyMesh::Point pp = mesh.data(*fit).center;
		double temp=(p-pp).length();
		if (temp<distance)
		{
			distance=temp;
			result=fit.handle();
		}
	}
	return result;
}

MyMesh::VertexHandle TriMesh::closestVertex(const Point &p)
{
	// find closest point;
	double distance=1.0e50;
	MyMesh::VertexHandle result;
	for(MyMesh::VertexIter vit = mesh.vertices_begin(); vit != mesh.vertices_end(); ++vit)
	{
		MyMesh::Point pp = mesh.point(vit);
		double temp=(p-pp).length();
		if (temp<distance)
		{
			distance=temp;
			result=vit.handle();
		}
	}
	return result;
}

double TriMesh::proc(const gmVector3 & x)
{
	return distancePointToMesh(Point((float)x[0],(float)x[1],(float)x[2]));
}

//TO DO: implement this: 
Intervald TriMesh::proc(const Box<double>& x)
{
	return Intervald();
}

double TriMesh::exactDistanceToMesh(const MyMesh::Point &p)
{
	double minDistance = 1000;
	double distance = 1000;
	MyMesh::FaceHandle fh;
	MyMesh::FaceHandle minFh;

	for(MyMesh::FaceIter fit = mesh.faces_begin(); fit != mesh.faces_end(); ++fit)
	{
		fh = fit.handle();
		
	    Point n=mesh.normal(fh);			
		distance = distancePointToTriangle(p, fh);
			
		if( fabs(distance) < fabs(minDistance) )
		{
				minDistance = distance;
				minFh = fh;
				
		}

	}
	
	MyMesh::FaceVertexIter fv_it=mesh.fv_iter(minFh);
	MyMesh::VertexHandle va=(fv_it.handle());
	
	Point pa=mesh.point(va);
	
	Point n=mesh.normal(minFh);
	gmVector3 pointVec(pa[0] - p[0], pa[1] - p[1], pa[2] - p[2]);
	gmVector3 normal(n[0], n[1], n[2]);

	if(dot(pointVec, normal) > 0)
		minDistance *= -1;
	
	return minDistance;
}

