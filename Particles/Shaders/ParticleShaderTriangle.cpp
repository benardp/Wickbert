#include "ParticleShaderTriangle.h"

REGISTER_PARTICLESTUFF(ParticleShaderTriangle,"Shader:ParticleShaderTriangle");

ParticleShaderTriangle::ParticleShaderTriangle(Particles *ps)
	:ParticleShader(ps,std::string("ParticleShaderTriangle"))
{
	new PSParamDouble(this,&lineWidth,1.0,"thickness","Line Thickness","If zero then draw faces, otherwise draw wireframe with desired thickness.");
	new PSParamBool(this,&faceNormal,true,"facenorm","Face Normal","true -> face normals, false -> vertex normals");

	new Attached<ParticleMesh>(this,&pm);
}

void ParticleShaderTriangle::drawPre()
{
	ParticleOpenMesh &mesh=pm->mesh;
	mesh.update_face_normals();
	if (!faceNormal)
		mesh.update_vertex_normals();
	// Draw triangle mesh
	if (lineWidth>0)
		drawWireframe(mesh);
	else
		drawPolygon(mesh);
}

void ParticleShaderTriangle::drawPolygon(ParticleOpenMesh &mesh)
{
	glPolygonOffset(1.0,1.0);
	glEnable(GL_POLYGON_OFFSET_FILL);
	glBegin(GL_TRIANGLES);
	for (ParticleOpenMesh::FaceIter f_it=mesh.faces_begin(); f_it!=mesh.faces_end();++f_it)
	{
		if (faceNormal) {
			//glNormal3fv(f_it->normal());
			//glNormal3fv(mesh.normal(f_it)); //changed by Matei Stroila
			ParticleOpenMesh::Normal n = mesh.normal(f_it);
			glNormal3f(n[0],n[1],n[2]);
		}
		for (ParticleOpenMesh::FaceVertexIter fv_it=mesh.fv_iter(f_it.handle());fv_it;++fv_it)
		{
			if (!faceNormal) {
				//glNormal3fv(fv_it->normal());
				//glNormal3fv(mesh.normal(fv_it)); //changed by Matei Stroila
				ParticleOpenMesh::Normal n = mesh.normal(fv_it);
				glNormal3f(n[0],n[1],n[2]);
			}
			//glVertex3fv(mesh.point(fv_it));//changed by Matei Stroila (*vit)
			ParticleOpenMesh::Point p = mesh.point(fv_it);
			glVertex3f(p[0],p[1],p[2]);
		}
	}
	glEnd();
	glDisable(GL_POLYGON_OFFSET_FILL);
}

void ParticleShaderTriangle::drawWireframe(ParticleOpenMesh &mesh)
{
	glDisable(GL_LIGHTING);
	for (ParticleOpenMesh::FaceIter f_it=mesh.faces_begin(); f_it!=mesh.faces_end();++f_it)
	{
		glBegin(GL_LINE_LOOP);
		for (ParticleOpenMesh::FaceVertexIter fv_it=mesh.fv_iter(f_it.handle());fv_it;++fv_it) {
			//glVertex3fv(mesh.point(fv_it));//changed by Matei Stroila (*vit)
			ParticleOpenMesh::Point p = mesh.point(fv_it);
			glVertex3f(p[0],p[1],p[2]);
		}
		glEnd();
	}
	glEnable(GL_LIGHTING);
}
