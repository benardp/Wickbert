#include "ParticleShaderCluster.h"
#include "SurfaceInterrogator.h"
#include "Surface/OpenMesh/ClusterMesh.h"

REGISTER_PARTICLESTUFF(ParticleShaderCluster,"Shader:ParticleShaderCluster");

ParticleShaderCluster::ParticleShaderCluster(Particles *ps)
	:ParticleShader(ps,std::string("ParticleShaderCluster"))
{
	lineWidth=1;
	useMaterial=1;
	drawOutline=0;
}

void ParticleShaderCluster::attachAttributes()
{
	// find the cluster mash
	attachAttribute(surInt,std::string("SurfaceInterrogator"));

	ParticleShader::attachAttributes();
}

int ParticleShaderCluster::qlen()
{
	return 3;
}

void ParticleShaderCluster::getq(double *q)
{
	q[0] = lineWidth;
	q[1] = useMaterial;
	q[2] = drawOutline;
}

void ParticleShaderCluster::setq(double *q)
{
	lineWidth = q[0];
	useMaterial = q[1];
	drawOutline = q[2];
}

void ParticleShaderCluster::qname(char **qn)
{
	qn[0] = "Line Width";
	qn[1] = "Use Material";
	qn[2] = "Draw Outline";
}

void ParticleShaderCluster::draw(int s)
{
	// the first color means not assigned to a cluster
	static float materialIndex[16][4]={{0,0,0,1},
		{0.2,1,1,1},{1,0.2,1,1},{1,1,0.2,1},{0.2,1,0.2,1},{0.2,0.2,1},
        {0.5,1,1,1},{1,0.5,1,1},{1,1,0.5,1},{0.5,1,0.5,1},{0.5,0.5,1,1},
		{0.8,1,1,1},{1,0.8,1,1},{1,1,0.8,1},{0.8,0,0.8,1},{0.8,0.8,1,1}};
	static GLfloat spec[4] = {1.0,1.0,1.0,1.0};
	static GLfloat amb[4] = {0.0,0.0,0.0,1.0};

	if (surInt==NULL)
		return;
	ClusterMesh *cm=dynamic_cast<ClusterMesh *>(surInt->surface);
	if (cm==NULL)
		return;
	ClusterOpenMesh &mesh=cm->mesh;

	glPolygonOffset(1.0,1.0);
	glEnable(GL_POLYGON_OFFSET_FILL);

	if (useMaterial)
	{
		glMaterialfv(GL_BACK,GL_SPECULAR,spec);
		glMaterialfv(GL_BACK,GL_AMBIENT,amb);
		glMaterialfv(GL_FRONT,GL_SPECULAR,spec);
		glMaterialfv(GL_FRONT,GL_AMBIENT,amb);
	}
	else
		glDisable(GL_LIGHTING);

	for (ClusterOpenMesh::FaceIter f_it=mesh.faces_begin(); f_it!=mesh.faces_end();++f_it)
	{
		// to identify each object 
		glLoadName(f_it.handle().idx());
		glBegin(GL_TRIANGLES);
		// material selection
		if (useMaterial)
		{
			glMaterialfv(GL_FRONT,GL_DIFFUSE,materialIndex[(f_it->cluster%15)+1]);
			glMaterialfv(GL_BACK,GL_DIFFUSE,materialIndex[(f_it->cluster%15)+1]);
		}
		else
			glColor4fv(materialIndex[(f_it->cluster%15)+1]);
		for (ClusterOpenMesh::FaceVertexIter fv_it=mesh.fv_iter(f_it.handle());fv_it;--fv_it)
		{
			//glNormal3fv(-fv_it->normal());
			glNormal3fv(-mesh.normal(fv_it));//changed by Matei Stroila
			glVertex3fv(mesh.point(fv_it));//changed by Matei Stroila (*vit)
		}
		glEnd();
	}

	if (!useMaterial)
		glEnable(GL_LIGHTING);

	if (drawOutline)
	{
		glDisable(GL_LIGHTING);
		glColor3f(1,1,1);
		glDisable(GL_POLYGON_OFFSET_FILL);
		for (ClusterOpenMesh::FaceIter f_it=mesh.faces_begin(); f_it!=mesh.faces_end();++f_it)
		{
			glBegin(GL_LINE_LOOP);
			for (ClusterOpenMesh::FaceVertexIter fv_it=mesh.fv_iter(f_it.handle());fv_it;++fv_it)
				glVertex3fv(mesh.point(fv_it));//changed by Matei Stroila (*vit)
			glEnd();
		}
		glEnable(GL_LIGHTING);
	}
}
