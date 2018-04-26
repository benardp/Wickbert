#include "ShaderMeshRenderer.h"
#include "MeshInterrogator.h"
#include "Surface/OpenMesh/SurfaceMesh.h"

REGISTER_PARTICLESTUFF(ShaderMeshRenderer,"Shader:ShaderMeshRenderer");

ShaderMeshRenderer::ShaderMeshRenderer(Particles *ps)
	:ParticleShader(ps,std::string("ShaderMeshRenderer"))
{
	new Attached<MeshInterrogator>(this,&_mi);
	new PSParamBool(this, &_lighting,true,"Light","Lighting", "Activate GL Lighting");
	new PSParamBool(this, &_wireframe,false,"Wireframe","Wireframe Mode", "Show Mesh in Wireframe mode");
	_mesh=0;
	_glWireframeDisplayList=0;
	_glDisplayList=0;
}

ShaderMeshRenderer::~ShaderMeshRenderer()
{
	glDeleteLists(1,_glWireframeDisplayList);
	glDeleteLists(1,_glDisplayList);
}


void ShaderMeshRenderer::drawPost()
{
	if (!_mi) return;
	
	SurfaceMesh *mesh = dynamic_cast<SurfaceMesh *>(_mi->getMesh());
	if (!mesh) return;

	if (mesh!=_mesh)
	{
		_mesh=mesh;
		updateDisplayLists();
	}

	if (_lighting)
	{
		glEnable(GL_LIGHTING);
	}
	else 
	{
		glDisable(GL_LIGHTING);
	}				

	if (_wireframe)
	{
		glCallList(_glWireframeDisplayList);
	}
	else
		glCallList(_glDisplayList);

}

void ShaderMeshRenderer::updateDisplayLists()
{

	glPushAttrib(GL_ALL_ATTRIB_BITS);
	{
		glDeleteLists(1,_glWireframeDisplayList);
		_glWireframeDisplayList= glGenLists(1);
		glNewList(_glWireframeDisplayList,GL_COMPILE);
		{
			unsigned int i = 0;
			glColor3f(0.8f,0.8f,0.8f);

			glBegin(GL_LINES);

			for (SurfaceMesh::ConstEdgeIter e_it = _mesh->edges_begin(); e_it!=_mesh->edges_end(); ++e_it) 
			{
				SurfaceMesh::Point	p;
				SurfaceMesh::Normal n;
				n = _mesh->normal(_mesh->from_vertex_handle(_mesh->halfedge_handle(e_it.handle(),0)));
				p = _mesh->point(_mesh->from_vertex_handle(_mesh->halfedge_handle(e_it.handle(),0)));
				glNormal3f(n[0],n[1],n[2]);
				glVertex3f(p[0],p[1],p[2]);

				n = _mesh->normal(_mesh->to_vertex_handle(_mesh->halfedge_handle(e_it.handle(),0)));
				p = _mesh->point(_mesh->to_vertex_handle(_mesh->halfedge_handle(e_it.handle(),0)));
				glNormal3f(n[0],n[1],n[2]);
				glVertex3f(p[0],p[1],p[2]);

				/* This was the original code...
				glNormal3fv(_mesh->normal(_mesh->from_vertex_handle(_mesh->halfedge_handle(e_it.handle(),0))));
				glVertex3fv(_mesh->point(_mesh->from_vertex_handle(_mesh->halfedge_handle(e_it.handle(),0))));
				glNormal3fv(_mesh->normal(_mesh->to_vertex_handle(_mesh->halfedge_handle(e_it.handle(),0))));
				glVertex3fv(_mesh->point(_mesh->to_vertex_handle(_mesh->halfedge_handle(e_it.handle(),0))));
				...*/
			}
			glEnd();
		}
		glEndList();


		glDeleteLists(1,_glDisplayList);
		_glDisplayList= glGenLists(1);
		glNewList(_glDisplayList,GL_COMPILE);
		{
#if 0
			for (SurfaceMesh::ConstFaceIter f_it = _mesh->faces_begin(); f_it!=_mesh->faces_end(); ++f_it) 
			{
				int n=0;					
				glBegin(GL_TRIANGLES);
				for (SurfaceMesh::ConstFaceVertexIter fv_it = _mesh->cfv_iter(f_it.handle()); fv_it; ++fv_it) 
				{
					++n;
					glNormal3fv(_mesh->normal(fv_it.handle()));
					glVertex3fv(_mesh->point(fv_it.handle()));
				}
				glEnd();
				//if this assert went of, the class Surface Mesh contains a non-triangle type.
				//all that needs to be done is change the above "#if 1" to "#if 0"				
				assert(n==3);
			}
#else
			for (SurfaceMesh::ConstFaceIter f_it = _mesh->faces_begin(); f_it!=_mesh->faces_end(); ++f_it) 
			{
				glBegin(GL_POLYGON);
				for (SurfaceMesh::ConstFaceVertexIter fv_it = _mesh->cfv_iter(f_it.handle()); fv_it; ++fv_it) 
				{
					SurfaceMesh::Point	p;
					SurfaceMesh::Normal n;

					n = _mesh->normal(fv_it.handle());
					p = _mesh->point(fv_it.handle());

					glNormal3f(n[0],n[1],n[2]);
					glVertex3f(p[0],p[1],p[2]);

					/* original code...
					glNormal3fv(_mesh->normal(fv_it.handle()));
					glVertex3fv(_mesh->point(fv_it.handle()));
					...*/
				}
				glEnd();
			}
#endif
		}
		glEndList();
	}
	glPopAttrib();
}