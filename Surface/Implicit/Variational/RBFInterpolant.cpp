//RBFInterpolant.cpp
//Loads a mesh and interpolates it using radial basis functions placed at each of the vertices
//Mike Mullan 

#include "RBFInterpolant.h"
// implicit should not depend on particles, the surface should have its points stored in its own format
// particles are a way of displaying surfaces. This is a design flaw.
// #include "ParticlesInclude.h"


//Read in an OBJ file
bool RBFMesh::readFile(const char *filename)
{
	
	if (filename==NULL)
		return false;
	mesh.clear();
	OpenMesh::IO::read_mesh(mesh, filename);
	mesh.update_face_normals();
	mesh.update_vertex_normals();
	return true;

}

//Place RBF centers at each of the vertices
void RBFInterpolant::placeCenters()
{
			
  moving = 0;
  morph = 1;
  for(MyMesh::VertexIter vit = rbfMesh.mesh.vertices_begin(); vit != rbfMesh.mesh.vertices_end(); ++vit)
  {
		
		//We give every center a normal vector, even if we are using dipoles.  This is a hack so that we can use the same code
		//to operate on dipoles and regular on surface centers
	    
	 
		MyMesh::Point p = rbfMesh.mesh.point(vit);
		
		MyMesh::Point n = rbfMesh.mesh.normal(vit);
		gmVector3 point(p[0], p[1], p[2]);
		gmVector3 normal(n[0], n[1], n[2]);

		if(dipoles)
		{
			normals.push_back(normal);
			normals.push_back(normal);
			centers.push_back(RBFPoint(point - epsilon * normal, -epsilon));
			centers.push_back(RBFPoint(point + epsilon * normal, epsilon));
		}
		else
		{
			centers.push_back(RBFPoint(point, 0));
			normals.push_back(normal);
			
		}
  }

  //Place positive valued centers off the surface
  if(!dipoles)
  {
	 centers.push_back(RBFPoint(gmVector3(-5,0,0), 5));
	 centers.push_back(RBFPoint(gmVector3(5,0,0), 5));
	 centers.push_back(RBFPoint(gmVector3(0,-5,0), 5));
	 centers.push_back(RBFPoint(gmVector3(0,5,0), 5));
	 centers.push_back(RBFPoint(gmVector3(0,0,-5), 5));
	 centers.push_back(RBFPoint(gmVector3(0,0,5), 5));
	
  }


  updateRBF();

}


unsigned int RBFInterpolant::qlen()
{
	
	return ( centers.size() * 4 + 4 );
}


void RBFInterpolant::getq(double* q)
{
	
	unsigned int i = 0;
    for(i = 0; i < centers.size(); i++)
	{
		q[i * 4] = centers[i].c[0];
		q[i * 4 + 1] = centers[i].c[1];
		q[i * 4 + 2] = centers[i].c[2];
		q[i * 4 + 3] = centers[i].w;
	}

	q[i * 4] = m_p[0];
	q[i * 4 + 1] = m_p[1];
	q[i * 4 + 2] = m_p[2];
	q[i * 4 + 3] = m_p[3];
	
}

void RBFInterpolant::_setq(double* q)
{
	unsigned int i = 0;
	for(i = 0; i < centers.size(); i++)
	{
		centers[i].c[0] = q[i * 4];
		centers[i].c[1] = q[i * 4 + 1];
	    centers[i].c[2] = q[i * 4 + 2];
		centers[i].w = q[i * 4 + 3];
	}

    m_p[0] = q[i * 4];
	m_p[1] = q[i * 4 + 1];
	m_p[2] = q[i * 4 + 2];
	m_p[3] = q[i * 4 + 3];

}

void RBFInterpolant::procq(const gmVector3 & x, double* dfdq)
{
	unsigned int i = 0;
	for(i = 0; i < centers.size(); i++)
	{
		gmVector3 vec = centers[i].c - x;
		double r = vec.length();
		vec = vec * 3 * r * centers[i].w;
		dfdq[i * 4] = vec[0];
		dfdq[i * 4 + 1] = vec[1];
		dfdq[i * 4 + 2] = vec[2];
		dfdq[i * 4 + 3] = r * r * r;
	}
	dfdq[i * 4] = 1;
	dfdq[i * 4 + 1] = x[0];
	dfdq[i * 4 + 2] = x[1];
	dfdq[i * 4 + 3] = x[2];
}


void RBFInterpolant::getqname(char** qn) 
{
 
 
  for(unsigned int i = 0; i < qlen(); i++)
  {

	qn[i] = " ";

		//qn[i] = new char[4];
		//qn[i][0] = '\0';
	 //   sprintf(s,"%d",i);
		//strcat(qn[i],s);
	

  }

}

 //Hessian used for normal vector propagation
gmMatrix3 RBFInterpolant::computeHessian(gmVector3 point)
{
	gmMatrix3 hessian;
	gmMatrix3 hessianCenter;

	hessianCenter = hessianCenter.identity();

	for(unsigned int i = 0; i < centers.size(); i++)
	{
	
		gmVector3 diffs;
		diffs[0] = point[0] - centers[i].c[0];
		diffs[1] = point[1] - centers[i].c[1];
		diffs[2] = point[2] - centers[i].c[2];

		double r = sqrt(diffs[0] * diffs[0] + diffs[1] * diffs[1] + diffs[2] * diffs[2]);
		
		for(int row = 0; row < 3; row++)
		{
			for(int col = 0; col < 3; col++)
			{
				hessianCenter[row][col] *= r;
				hessianCenter[row][col] += diffs[row] * diffs[col] / r;
			}
		}
		
		hessianCenter *= 3 * centers[i].w;
        hessian += hessianCenter;
		
		hessianCenter = hessianCenter.identity();
	}
	
	
	return hessian;
	
}

// should moved to particle
//void RBFInterpolant::interpolate(Particles *ps, std::valarray<bool> flexible)
//{
//	SurfaceDeformation *sd = ps->getBehavior<SurfaceDeformation>(std::string("SurfaceDeformation"));
//    assert(sd);
//
//	ParticlePosition *position=ps->getAttribute<ParticlePosition>(std::string("ParticlePosition"));
//	assert(position);
//
//	ParticleVelocity *velocity=ps->getAttribute<ParticleVelocity>(std::string("ParticleVelocity"));
//	assert(velocity);
//	
//	//If the surface is moving, set the positions of the control particles to the positions of the RBF centers
//	if(moving)
//	{
//				
//		for(unsigned int i = 0; i < ps->size(); i ++)
//		{
//    		position->setPosition(i, centers[i].c);
//		}
//	
//		//adaptivity
//		if(ps->size() < centers.size())
//		{
//			int diff = centers.size() - ps->size();
//			if(dipoles)
//			{
//				ps->addParticle();
//				position->setPosition(ps->size() - 1, centers[normals.size() - 1].c);
//				ps->addParticle();
//				position->setPosition(ps->size() - 1, centers[normals.size() - 2].c);
//			}
//			else
//			{
//				
//				ps->addParticle();
//				gmVector3 temp = position->x[normals.size() - 1];
//                position->setPosition(ps->size() - 1, centers[normals.size() - 1].c);
//				position->x.back() = temp;
//			}
//
//
//			for(i = 0; i < diff; i++)
//			{
//				
//
//			}
//
//		}
//	}
//	//If the surface is stationary, set the positions of the RBF centers to the positions of the control particles
//	else
//	{
//		
//	//	centers.clear();
//	//	centers.resize(ps->size() + 1);
//		for (unsigned int i = 0; i < ps->size(); i++)
//	//		centers[i] = RBFPoint(position->getPosition(i),sd->h[i]);
//			centers[i].c = position->x[i];
//		update();
//	}
//
//	moving = 0;
//}

//Distance function to target mesh
double RBFInterpolant::targetShape(gmVector3 point)
{
	return targetMesh.exactDistanceToMesh(MyMesh::Point((float)point[0],(float)point[1],(float)point[2]));
}
