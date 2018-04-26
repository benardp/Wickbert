/** @file SurfaceMesh.h
	@author John C. Hart (based on TriMesh by Wen Yu Su
	@date 3 Jan. 2005
	SurfaceMesh is a kind of Surface that contains an OpenMesh mesh.
*/

#ifndef SURFACEMESH_H
#define SURFACEMESH_H

//wxGTK got conflicting types in OpenMesh, so undefed  them here -ms
#ifdef __WXGTK__
#ifdef True 
#undef True
#endif
#ifdef False 
#undef False
#endif
#ifdef None
#undef None
#endif
#ifdef Status
#undef Status
#endif
#endif

#ifdef check
#undef check
#endif

// openmesh includes
#include "OpenMesh/Core/IO/BinaryHelper.hh"
#include "OpenMesh/Core/IO/MeshIO.hh"
#include "OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh"
#ifdef _MSC_VER
#  ifndef OM_STATIC_BUILD
#  define OM_STATIC_BUILD
#  endif
#  define INCLUDE_TEMPLATES
#  include <OpenMesh/Core/IO/IOInstances.hh>
#endif

#include <vector>
#include <string>

#include "Surface/Surface.h"
#include "libgm/gmVec3.h"


using namespace std;
using namespace OpenMesh;

struct SurfaceMeshTraits : public OpenMesh::DefaultTraits
{
	FaceAttributes(OpenMesh::Attributes::Normal);
	VertexAttributes(OpenMesh::Attributes::Normal | OpenMesh::Attributes::TexCoord2D);
};

typedef OpenMesh::PolyMesh_ArrayKernelT<SurfaceMeshTraits> MyMesh;

class SurfaceMesh : public Surface, public OpenMesh::PolyMesh_ArrayKernelT<SurfaceMeshTraits>
{
public:
	SurfaceMesh(std::string fname = "") : filename(fname)
	{
		if ( readFile(filename) )
			objectName = filename;
		else
			objectName = std::string("Not found: ") + filename;

		registry_name = "OpenMesh:PolyMesh:ArrayKernel";
	}


	bool readFile(std::string fname)
	{
		clear();
		if ( fname == "" ) return false;
		filename = fname;
		if ( OpenMesh::IO::read_mesh(*this, fname) ) {
			scaleMesh();
			update_face_normals();
			update_vertex_normals();
			return true;
		} else
			return false;
	}


	bool readStream(std::istream& is)
	{
		string mname, p, fname;
		char c;
		bool status(true);
		
		is >> mname;
		setObjectName(mname);
		is >> p;	// "{"

		is >> c;
		if (c != '"') {
			std::cerr << "Specified mesh " << mname << "does not have a quoted filename.\n";
			status = false;
		}

		is >> std::noskipws;
		while ((is >> c) && c != '"') fname += c;
		is >> std::skipws;

		if ( ! readFile(fname) )
			status = false;

		registry_name = "OpenMesh:PolyMesh:ArrayKernel";

		is >> p;	// "}"

		return status;
	}


	bool writeStream(std::ostream& out)
	{
		out << objectName << " {\n";
		out << "\t\"" << filename << "\"\n";
		out << "}\n\n";

		return true;
	}


	MyMesh::Point bbox_min, bbox_max;
	
	void computeBoundingBox()
	{
		MyMesh::VertexIter vit = this->vertices_begin(), vend = this->vertices_end();

		MyMesh::Point firstPoint = this->point(vit);
		bbox_min = bbox_max = firstPoint;

		for ( ; vit != vend; ++vit)
		{		
			MyMesh::Point p = this->point(vit);
			for (int i=0; i < 3; i++)
			{
				bbox_min[i] = ((bbox_min[i]<p[i])?bbox_min[i]:p[i]);
				bbox_max[i] = ((bbox_max[i]>p[i])?bbox_max[i]:p[i]);
			}
		}
	}


	void scaleMesh() {

		computeBoundingBox();

		MyMesh::Point bbox_mid;

		for (int i=0; i < 3; i++)
			bbox_mid[i] = 0.5f * (bbox_max[i] + bbox_min[i]);

		MyMesh::Point t = bbox_max - bbox_min;
		MyMesh::Scalar radius = 2.0f / (0.5f * sqrt(t[0] * t[0] + t[1] * t[1] + t[2] * t[2])); 
		MyMesh::VertexIter vit = this->vertices_begin(), vend = this->vertices_end();
		MyMesh::Point foo;

		for ( ; vit != vend; ++vit)
		{	
			foo = this->point(vit);

			foo = radius * (foo - bbox_mid);
			this->set_point(vit, foo);
		
		}
	}

	bool saveAs(std::string filename) {
		return OpenMesh::IO::write_mesh(*this,filename,IO::Options::VertexTexCoord);
	}

protected:
	std::string filename;
};

#endif
