/**
 * Implementation of a volume data implicit 
 * @file ITKImplicit.h
 * @date 10/01/2005
 * @author Matei N. Stroila modified a lot by Elmar Eisemann :)
 */

#ifdef WB_USE_ITK
#include "ITKImplicit.h"
#include "itkPoint.h"
#include "Particles/Attributes/ParticleBoundingBox.h"
#include "Surface/OpenMesh/SurfaceMesh.h"
//to calculate distance fields
#include "Surface/Implicit/DistanceField/signeddistancetransform.h"
#include <fstream>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
typedef itk::Point<double, Dimension> PointType;

REGISTER_IMPLICIT(ITKImplicit,"ITKImplicit");

ITKImplicit::ITKImplicit(void)
{
	new SurfParamString(this,&dirFile,"","fileDir", "Input file or dicom directory");
	new SurfParamButton(this,new LoadButtonCallback(this),"load","","load a dicom or imagetype file (extension .vtk suggested)");
	new SurfParamButton(this,new SaveButtonCallback(this),"save","","save a dicom or imagetype file (extension .vtk suggested)");
	new SurfParamDouble(this,&threshold,0,"Threshold", "CT Number");
	new SurfParamDouble(this,&scale,1.0,"Scale", "Volume Scaling Factor");
	new SurfParamgmVector3(this,&tr,gmVector3(0,0,0),"Translation","Volume Translation");
	new SurfParamButton(this, new ApplyParameterChangesCallback(this),"Apply","Apply parameter changes",  "apply parameter changes: Namely translation, scaling,"
																											"and spline order");

	new SurfSurfRefParam(this, &surface,"<invalid/optional>","Surf(opt.)","SurfaceReference(optional)",	
																			"The ITK Volume can also be initialized using a" 
																			"surface. In the case of an implicit the bounding box should be meaningful."
																			"If no bounding box is provided a unit cube centered at the origin is assumed."
																			"The implicit is sampled inside of this bounding box."
																			"In the case of a surface mesh, the behavior is special:"
																			"a distance field is calculated inside of ITS bounding box");


	new SurfParamInt(this,&(voxelRes[0]),64,"sizeX", "NbX", "Nb of gridcells along x-Axis");
	new SurfParamInt(this,&(voxelRes[1]),64,"sizeY","nbY", "Nb of gridcells along y-Axis");
	new SurfParamInt(this,&(voxelRes[2]),64,"sizeZ","nbZ", "Nb of gridcells along z-Axis");
	new SurfParamDouble(this, &boundingBoxScaling, 1.3,"scale","scalebb", "Scaling of the bounding box to make it slightly bigger");

	new SurfAttrRefParam(this,(ParticleAttribute **)&bbox,"invalid:ParticleBoundingBox","bbox","box",
		"Particle bounding box attribute (it is advised to provide one for implicit surfaces).");

	new SurfParamButton(this, new SetITKImplicit(this),"initImp","Initialize ITK implicit","set ITK imp");
	new SurfParamInt(this,&order,3,"Order","Spline Order");


	//DICOM reader
	dicomIO = ImageIOType::New();
	reader = ReaderType::New();
	reader->SetImageIO( dicomIO );
	nameGenerator = NamesGeneratorType::New();


	//now we initialize a simple default image
	//and all the rest of this constructor is just to do this...
	//...I love the ITK interface ... lol
	//init default image, this can still be overwritten
	myImage = ImageType::New();

	ImageType::SpacingType spacing;
	spacing[0]=0.1;
	spacing[1]=0.1;;
	spacing[2]=0.1;;
	myImage->SetSpacing(spacing);
	
	//we set the origin such that center lies at (0,0,0)
	ImageType::PointType origin;
	origin[0]=-0.05;
	origin[1]=-0.05;
	origin[2]=-0.05;
	myImage->SetOrigin(origin);
	
	//we want an image with voxelRes voxels
	ImageType::IndexType start;
	start[0]=0;
	start[1]=0;
	start[2]=0;
	ImageType::SizeType size;
	size[0]=1;
	size[1]=1;
	size[2]=1;
	ImageType::RegionType region;
	region.SetIndex(start);
	region.SetSize(size);
	myImage->SetRegions(region);
	
	//allocate memory
	myImage->Allocate();


	//interpolating function
	myInterpFunc = InterpFunc::New();	
	myInterpFunc->SetSplineOrder((unsigned int) 0);
	myInterpFunc->SetInputImage(myImage);

	
}

double ITKImplicit::proc(const gmVector3 & x)
{	
	ImageType::PointType myPoint;
	myPoint[0]=x[0];
	myPoint[1]=x[1];
	myPoint[2]=x[2];

	if(myInterpFunc->IsInsideBuffer(myPoint))	
		return myInterpFunc->Evaluate(myPoint) - threshold;
	else
	{
		//we want to attract particles in direction of the center
		//later we might want to add that particles slow down before entering the volume
		//for the moment we just use this simple formula
		ImageType::SpacingType spacing = myImage->GetSpacing();
		ImageType::SizeType size = myImage->GetLargestPossibleRegion().GetSize();


		gmVector3 scaling;
		scaling[0] = 1.0/(spacing[0]*size[0]/2.0);
		scaling[1] = 1.0/(spacing[1]*size[1]/2.0);
		scaling[2] = 1.0/(spacing[2]*size[2]/2.0);

		scaling=scaling.normalize();

		gmVector3 ellipse=(x-center);

		ellipse[0]*=scaling[0];
		ellipse[1]*=scaling[1];
		ellipse[2]*=scaling[2];

		return 1;//ellipse.lengthSquared();
	}
}


gmVector3 ITKImplicit::grad(const gmVector3 & x)
{
	ImageType::PointType myPoint;
	myPoint[0]=x[0];
	myPoint[1]=x[1];
	myPoint[2]=x[2];

	if(myInterpFunc && myInterpFunc->IsInsideBuffer(myPoint))
	{	
		itk::CovariantVector<double> v(myInterpFunc->EvaluateDerivative(myPoint));
		//before it was dependent on the scaling, but I see no reason for it...
		return gmVector3(v[0], v[1], v[2]);
	}
	else
	{
		//we return a grad which is like one on the sphere to attract particles to the center...
		ImageType::SpacingType spacing = myImage->GetSpacing();
		ImageType::SizeType size = myImage->GetLargestPossibleRegion().GetSize();


		gmVector3 scaling;
		scaling[0] = 1.0/(spacing[0]*size[0]/2.0);
		scaling[1] = 1.0/(spacing[1]*size[1]/2.0);
		scaling[2] = 1.0/(spacing[2]*size[2]/2.0);
		scaling=scaling.normalize();



		gmVector3 ellipse=(x-center);

		ellipse[0]*=scaling[0];
		ellipse[1]*=scaling[1];
		ellipse[2]*=scaling[2];
		
		return (ellipse).normalize(); 
	}
}


gmMatrix3 ITKImplicit::hess(const gmVector3 & x)
{

	ImageType::PointType myPoint;
	myPoint[0]=x[0];
	myPoint[1]=x[1];
	myPoint[2]=x[2];

	if(myInterpFunc && myInterpFunc->IsInsideBuffer(myPoint))
	{
		itk::CovariantVector<double> grad(myInterpFunc->EvaluateDerivative(myPoint));

		ImageType::PointType myPointxEps(myPoint);
		myPointxEps[0]+=m_epsilon;
		itk::CovariantVector<double> gradxEps;
		if (myInterpFunc->IsInsideBuffer(myPointxEps))
			gradxEps=(myInterpFunc->EvaluateDerivative(myPointxEps));
		else
			return gmMatrix3(1,0,0,
							0,1,0,
							0,0,1); 

		gmVector3 gradx=(gmVector3(gradxEps[0], gradxEps[1],gradxEps[2])-gmVector3(grad[0],grad[1],grad[2]))/m_epsilon;


		ImageType::PointType myPointyEps(myPoint);
		myPointyEps[1]+=m_epsilon;
		itk::CovariantVector<double> gradyEps;
		if (myInterpFunc->IsInsideBuffer(myPointyEps))
			gradyEps=(myInterpFunc->EvaluateDerivative(myPointyEps));
		else
			return gmMatrix3(1,0,0,
							0,1,0,
							0,0,1); 

		gmVector3 grady=(gmVector3(gradyEps[0], gradyEps[1],gradyEps[2])-gmVector3(grad[0],grad[1],grad[2]))/m_epsilon;

		ImageType::PointType myPointzEps(myPoint);
		myPointzEps[2]+=m_epsilon;
		itk::CovariantVector<double> gradzEps;
		if (myInterpFunc->IsInsideBuffer(myPointyEps))
			gradzEps=(myInterpFunc->EvaluateDerivative(myPointzEps));
		else
			return gmMatrix3(1,0,0,
							0,1,0,
							0,0,1); 

		gmVector3 gradz=(gmVector3(gradzEps[0], gradzEps[1],gradzEps[2])-gmVector3(grad[0],grad[1],grad[2]))/m_epsilon;


		const double	fxx=gradx[0],	fxy=(grady[0]+gradx[1])/2.0,	fxz=(gradz[0]+gradx[2])/2.0,
										fyy=grady[1],					fyz=(grady[2]+gradz[1])/2.0,
																		fzz=gradz[2];

		return gmMatrix3(fxx,fxy,fxz,fxy,fyy,fyz,fxz,fyz,fzz);
	}
	else
	{
		//we return a hessian which is like one on the sphere to attract particles to the center...
		return gmMatrix3(1,0,0,
						 0,1,0,
						 0,0,1); 
	}
}




void ITKImplicit::load()
{
	bool success=true;
	if(dirFile.size() !=0)
	{
		//Apparently ther is no standard C++ way to test for a file existence.
		//It seems to me, that I rememeber something like this for system programming
		//using stats... but it has been a while and this should work too.
		std::ifstream file(dirFile.c_str());
		bool noFile=!file.is_open();
		file.close();

		if (noFile)
		{
			success=loadDicom(dirFile);
		}
		else
		{
			//it has to be an itk image file.
			//because we do not support anything else for the moment...
			itk::ImageFileReader<ImageType>::Pointer reader=itk::ImageFileReader<ImageType>::New();
				reader->SetFileName(dirFile.c_str());
			
			myImage=reader->GetOutput();
			
			try
			{
				reader->Update();
			}
			catch(itk::ExceptionObject &ex)
			{
				std::cerr << "Exception thrown while reading the image" << std::endl;
				std::cerr << ex << std::endl;

				success=false;
			}
		}
	}

	if (success)
	{

		originalSpacing = myImage->GetSpacing();
		originalOrigin = myImage->GetOrigin();
		//we apply eventual transformations 
		//and init the spline interpolator
		//we refit the particles bounding box
		applyParameterChanges();
	}
}
void ITKImplicit::save() const
{
	if (dirFile.size()==0)
		return;
	itk::ImageFileWriter<ImageType>::Pointer writer=itk::ImageFileWriter<ImageType>::New();
	writer->SetInput(myImage);
	writer->SetFileName(dirFile.c_str());
	try{
		writer->Update();
	}
	catch(itk::ExceptionObject &ex)
	{
		std::cerr << "Exception thrown while writing the image" << std::endl;
		std::cerr << ex << std::endl;
	}
}

void ITKImplicit::init(void)
{
	
	//init from surface
	//Via dynamic cast we test whether the surface type allows for a transformation
	//into an ITKImplicit
	if (!surface)
		return;
	
	//init default image, this can still be overwritten
	myImage = ImageType::New();
	
	SurfaceMesh * mesh=dynamic_cast<SurfaceMesh*>(surface);
	if (mesh)
	{
		mesh->computeBoundingBox();
		gmVector3 maxV(mesh->bbox_max[0], mesh->bbox_max[1], mesh->bbox_max[2]);
		gmVector3 minV(mesh->bbox_min[0], mesh->bbox_min[1], mesh->bbox_min[2]);
		initScalingAndOriginFromBoundingBox(minV,maxV);

		//we have a special surface, a mesh surface
		initFromSurfaceMesh(mesh);
	}
	//surface could be an implicit
	Implicit * imp=dynamic_cast<Implicit*>(surface);
	if (imp)
	{
		//if we have an implicit we want to check whether the particle bounding box is set.
		//if so, we adjust the image 
		if (bbox)
		{
			//we have a bounding box, so we use it to evaluate the
			//values.
			initScalingAndOriginFromBoundingBox(bbox->min, bbox->max);
		}
		else
		{
			initScalingAndOriginFromBoundingBox(gmVector3(-0.5,-0.5,-0.5),
												gmVector3(0.5,0.5,0.5));
		}
		initFromImplicit(imp);
	}

	originalSpacing = myImage->GetSpacing();
	originalOrigin = myImage->GetOrigin();

	//we apply eventual transformations 
	//and init the spline interpolator
	//we refit the particles bounding box
	applyParameterChanges();
}


void ITKImplicit::initScalingAndOriginFromBoundingBox(const gmVector3 & minV, const gmVector3 & maxV)
{
	//we make the bounding box slightly bigger
	
	gmVector3 originalCenter=(minV+maxV)/2;
	gmVector3 newMinV=(minV-originalCenter)*boundingBoxScaling+originalCenter;
	gmVector3 newMaxV=(maxV-originalCenter)*boundingBoxScaling+originalCenter;

	//for the moment we will use a scaling, 
	//that leads to a unit cube
	ImageType::SpacingType spacing;
	gmVector3 lengthV=newMaxV-newMinV;
	spacing[0]=fabs(lengthV[0])/((double)voxelRes[0]);
	spacing[1]=fabs(lengthV[1])/((double)voxelRes[1]);
	spacing[2]=fabs(lengthV[2])/((double)voxelRes[2]);
	myImage->SetSpacing(spacing);
	
	//we set the origin such that center lies at (0,0,0)
	ImageType::PointType origin;
	origin[0]=newMinV[0];
	origin[1]=newMinV[1];
	origin[2]=newMinV[2];
	myImage->SetOrigin(origin);
	
	//we want an image with voxelRes voxels
	ImageType::IndexType start;
	start[0]=0;
	start[1]=0;
	start[2]=0;
	ImageType::SizeType size;
	size[0]=max(1,voxelRes[0]);
	size[1]=max(1,voxelRes[1]);
	size[2]=max(1,voxelRes[2]);
	ImageType::RegionType region;
	region.SetIndex(start);
	region.SetSize(size);
	myImage->SetRegions(region);
	
	//allocate memory
	myImage->Allocate();
}

void ITKImplicit::applyParameterChanges()
{
	myInterpFunc->SetSplineOrder((unsigned int) order);
	myInterpFunc->SetInputImage(myImage);

	ImageType::SpacingType spacing;
	spacing[0]=originalSpacing[0]*scale;
	spacing[1]=originalSpacing[1]*scale;
	spacing[2]=originalSpacing[2]*scale;

	myImage->SetSpacing(spacing);


	ImageType::PointType origin;
	origin[0]=originalOrigin[0]+tr[0];
	origin[1]=originalOrigin[1]+tr[1];
	origin[2]=originalOrigin[2]+tr[2];

	myImage->SetOrigin(origin);

	
	
	ImageType::SizeType size;
	size = myImage->GetLargestPossibleRegion().GetSize();

	center[0]=origin[0]+0.5*spacing[0]*size[0];
	center[1]=origin[1]+0.5*spacing[1]*size[1];
	center[2]=origin[2]+0.5*spacing[2]*size[2];

	//finally we adjust the bounding box of the particles
	if(!bbox) return;

	bbox->min[0] = origin[0];
	bbox->min[1] = origin[1];
	bbox->min[2] = origin[2];
	bbox->max[0] = origin[0]+size[0]*spacing[0];
	bbox->max[1] = origin[1]+size[1]*spacing[1];
	bbox->max[2] = origin[2]+size[2]*spacing[2];
}



bool ITKImplicit::loadDicom(const std::string & directory)
{
	nameGenerator->SetInputDirectory( directory.c_str() );
	
	typedef std::vector< std::string >   FileNamesContainer;
	FileNamesContainer fileNames;
	fileNames = nameGenerator->GetInputFileNames( );
	unsigned int numberOfFilenames =  fileNames.size();

	for(unsigned int fni = 0; fni<numberOfFilenames; fni++)
	{
		std::cout << "filename # " << fni << " = ";
		std::cout << fileNames[fni] << std::endl;
	}
	
	reader->SetFileNames( fileNames );
		
	try
	{
		reader->Update();		
	}
	catch (itk::ExceptionObject &excp)
	{
		std::cerr << "Exception thrown while reading the image" << std::endl;
		std::cerr << excp << std::endl;
	
		return false;
	}
	
	// 
	// At this point we would have the volumetric data loaded in memory and we can
	// get access to it by invoking the \code{GetOutput()} method in the reader.
	//
	myImage = reader->GetOutput();		
	return true;
}

void ITKImplicit::initFromImplicit(Implicit * imp)
{
	ImageType::SizeType size = myImage->GetLargestPossibleRegion().GetSize();
	ImageType::PointType origin = myImage->GetOrigin();
	ImageType::SpacingType spacing = myImage->GetSpacing(); 

	ImageType::IndexType index;
	for (unsigned int z=0;z<size[2];++z)
	{
		index[2]=z;
		for (unsigned int y=0;y<size[1];++y)
		{
			index[1]=y;
			for (unsigned int x=0;x<size[0];++x)
			{
				//transform indices in index
				index[0]=x;

				//calculate world position of voxel
				gmVector3 point(origin[0]+x*spacing[0], origin[1]+y*spacing[1], origin[2]+z*spacing[2]);
				double value = imp->proc(point);

				//store value in image
				myImage->SetPixel(index, value);
			}
		}
	}
}







void ITKImplicit::initFromSurfaceMesh(SurfaceMesh * mesh)
{
	SignedDistanceTransform dt(*mesh);

	//fill up voxel grid values
	ImageType::SizeType size = myImage->GetLargestPossibleRegion().GetSize();
	ImageType::PointType origin = myImage->GetOrigin();
	ImageType::SpacingType spacing = myImage->GetSpacing(); 
	ImageType::IndexType index;

	for (unsigned int z=0;z<size[2];++z)
	{
		index[2]=z;
		for (unsigned int y=0;y<size[1];++y)
		{
			index[1]=y;
			for (unsigned int x=0;x<size[0];++x)
			{
				index[0]=x;
				myImage->SetPixel(index,DBL_MAX);
			}
		}
	}

	for (unsigned int z=0;z<size[2];++z)
	{
		index[2]=z;
		for (unsigned int y=0;y<size[1];++y)
		{
			index[1]=y;
			for (unsigned int x=0;x<size[0];++x)
			{
				index[0]=x;
				
				//calculate world position of voxel
				gmVector3 worldPoint(origin[0]+x*spacing[0], origin[1]+y*spacing[1], origin[2]+z*spacing[2]);
				myImage->SetPixel(index,dt.signedDistanceTo(worldPoint));
			}
		}
	}


	//Implementation using only the nearest point class
	//std::vector<TriangleNearestPoint> distanceTester;
	//std::vector<gmVector3> normals;

	//for (SurfaceMesh::ConstFaceIter f_it = mesh->faces_begin(); f_it!=mesh->faces_end(); ++f_it) 
	//{
	//	SurfaceMesh::Point n(mesh->normal(f_it.handle()));
	//	normals.push_back(gmVector3(n[0],n[1],n[2]));

	//	std::vector<gmVector3> vertices;
	//	for (SurfaceMesh::ConstFaceVertexIter fv_it = mesh->cfv_iter(f_it.handle()); fv_it; ++fv_it) 
	//	{
	//		SurfaceMesh::Point p= mesh->point(fv_it.handle());
	//		vertices.push_back(gmVector3(p[0],p[1],p[2]));

	//	}
	//	//test if mesh is triangular
	//	if (vertices.size()>3)
	//		return;
	//	assert(vertices.size()==3);
	//	distanceTester.push_back(TriangleNearestPoint(vertices[0],vertices[1],vertices[2]));
	//}

	////fill up voxel grid values
	//ImageType::SizeType size = myImage->GetLargestPossibleRegion().GetSize();
	//ImageType::PointType origin = myImage->GetOrigin();
	//ImageType::SpacingType spacing = myImage->GetSpacing(); 

	//TriangleNearestPoint::VoronoiRegion region;
	//gmVector3 nearestPoint;
	//ImageType::IndexType index;
	////init image values to FLT_MAX
	//for (unsigned int z=0;z<size[2];++z)
	//{
	//	index[2]=z;
	//	for (unsigned int y=0;y<size[1];++y)
	//	{
	//		index[1]=y;
	//		for (unsigned int x=0;x<size[0];++x)
	//		{
	//			index[0]=x;
	//			myImage->SetPixel(index,DBL_MAX);
	//		}
	//	}
	//}

	//for (unsigned int z=0;z<size[2];++z)
	//{
	//	index[2]=z;
	//	for (unsigned int y=0;y<size[1];++y)
	//	{
	//		index[1]=y;
	//		for (unsigned int x=0;x<size[0];++x)
	//		{
	//			index[0]=x;
	//			
	//			//calculate world position of voxel
	//			gmVector3 worldPoint(origin[0]+x*spacing[0], origin[1]+y*spacing[1], origin[2]+z*spacing[2]);
	//			
	//			for (unsigned int i=0;i<distanceTester.size();++i)
	//			{
	//				nearestPoint=distanceTester[i].NearestPointTo(worldPoint,region);
	//				gmVector3 diff=(worldPoint-nearestPoint);
	//				double value=diff.length();
	//				double oldValue=myImage->GetPixel(index);
	//				if (fabs(oldValue)>value)
	//				{
	//					value=(dot(diff,normals[i])>0)?value:-value;
	//					myImage->SetPixel(index,value);	
	//				}
	//			}
	//		}
	//	}
	//}
}



#endif




//DAVE EBERLY'S version for the distance calculation
//In our system...
//Not tested yet!
#if 0
	TriangleNearestPoint::VoronoiRegion region;
	gmVector3 nearestPoint;
	ImageType::IndexType index;
	//init image values to FLT_MAX
	for (unsigned int z=0;z<size[2];++z)
	{
		index[2]=z;
		for (unsigned int y=0;y<size[1];++y)
		{
			index[1]=y;
			for (unsigned int x=0;x<size[0];++x)
			{
				index[0]=x;
				myImage->SetPixel(index,DBL_MAX);
			}
		}
	}


	for (unsigned int z=0;z<size[2];++z)
	{
		index[2]=z;
		for (unsigned int y=0;y<size[1];++y)
		{
			index[1]=y;
			for (unsigned int x=0;x<size[0];++x)
			{
				index[0]=x;
				
				//calculate world position of voxel
				gmVector3 worldPoint(origin[0]+x*spacing[0], origin[1]+y*spacing[1], origin[2]+z*spacing[2]);
				for (SurfaceMesh::ConstFaceIter f_it = mesh->faces_begin(); f_it!=mesh->faces_end(); ++f_it) 
				{
					std::vector<gmVector3> vertices;
					for (SurfaceMesh::ConstFaceVertexIter fv_it = mesh->cfv_iter(f_it.handle()); fv_it; ++fv_it) 
					{
						vertices.push_back(gmVector3(p[0],p[1],p[2]));
					}
					gmVector3 nearest;
					double dist=fabs(distancePointTriangle(vertices[0],vertices[1],vertices[2],worldPoint, &nearest));
					double oldDist=myImage->GetPixel(index);
					if (fabs(dist)<fabs(oldDist))
					{
						SurfaceMesh::Point normal=mesh->normal(f_it.handle);
						gmVector3 n(normal[0],normal[1],normal[2]);
						dist=(dot(nearest,(worldPoint-vertices[0]))>0
						myImage->SetPixel(index,dist);
					}
				}
			}
		}
	}


	
	
double distancePointTriangle(const gmVector3 & v0, const gmVector3 & v1, const gmVector3 & v2,const gmVector3 &p, gmVector3 * nearest)
{
	assert(nearest);

    gmVector3 kDiff = v0 - p;
    gmVector3 kEdge0 = v1 - v0;
    gmVector3 kEdge1 = v2 - v0;
    double fA00 = kEdge0.lengthSquared();
    double fA01 = dot(kEdge0,kEdge1);
    double fA11 = kEdge1.lengthSquared();
    double fB0 = dot(kDiff,kEdge0);
    double fB1 = dot(kDiff,kEdge1);
    double fC = kDiff.lengthSquared();
    double fDet = fabs(fA00*fA11-fA01*fA01);
    double fS = fA01*fB1-fA11*fB0;
    double fT = fA01*fB0-fA00*fB1;
    double fSqrDistance;

    if (fS + fT <= fDet)
    {
        if (fS < 0.0)
        {
            if (fT < 0.0)  // region 4
            {
                if (fB0 < 0.0)
                {
                    fT = 0.0;
                    if (-fB0 >= fA00)
                    {
                        fS = 1.0;
                        fSqrDistance = fA00+(2.0)*fB0+fC;
                    }
                    else
                    {
                        fS = -fB0/fA00;
                        fSqrDistance = fB0*fS+fC;
                    }
                }
                else
                {
                    fS = 0.0;
                    if (fB1 >= 0.0)
                    {
                        fT = 0.0;
                        fSqrDistance = fC;
                    }
                    else if (-fB1 >= fA11)
                    {
                        fT = 1.0;
                        fSqrDistance = fA11+(2.0)*fB1+fC;
                    }
                    else
                    {
                        fT = -fB1/fA11;
                        fSqrDistance = fB1*fT+fC;
                    }
                }
            }
            else  // region 3
            {
                fS = 0.0;
                if (fB1 >= 0.0)
                {
                    fT = 0.0;
                    fSqrDistance = fC;
                }
                else if (-fB1 >= fA11)
                {
                    fT = 1.0;
                    fSqrDistance = fA11+(2.0)*fB1+fC;
                }
                else
                {
                    fT = -fB1/fA11;
                    fSqrDistance = fB1*fT+fC;
                }
            }
        }
        else if (fT < 0.0)  // region 5
        {
            fT = 0.0;
            if (fB0 >= 0.0)
            {
                fS = 0.0;
                fSqrDistance = fC;
            }
            else if (-fB0 >= fA00)
            {
                fS = 1.0;
                fSqrDistance = fA00+(2.0)*fB0+fC;
            }
            else
            {
                fS = -fB0/fA00;
                fSqrDistance = fB0*fS+fC;
            }
        }
        else  // region 0
        {
            // minimum at interior point
            double fInvDet = (1.0)/fDet;
            fS *= fInvDet;
            fT *= fInvDet;
            fSqrDistance = fS*(fA00*fS+fA01*fT+(2.0)*fB0) +
                fT*(fA01*fS+fA11*fT+(2.0)*fB1)+fC;
        }
    }
    else
    {
        double fTmp0, fTmp1, fNumer, fDenom;

        if (fS < 0.0)  // region 2
        {
            fTmp0 = fA01 + fB0;
            fTmp1 = fA11 + fB1;
            if (fTmp1 > fTmp0)
            {
                fNumer = fTmp1 - fTmp0;
                fDenom = fA00-2.0f*fA01+fA11;
                if (fNumer >= fDenom)
                {
                    fS = 1.0;
                    fT = 0.0;
                    fSqrDistance = fA00+(2.0)*fB0+fC;
                }
                else
                {
                    fS = fNumer/fDenom;
                    fT = 1.0 - fS;
                    fSqrDistance = fS*(fA00*fS+fA01*fT+2.0f*fB0) +
                        fT*(fA01*fS+fA11*fT+(2.0)*fB1)+fC;
                }
            }
            else
            {
                fS = 0.0;
                if (fTmp1 <= 0.0)
                {
                    fT = 1.0;
                    fSqrDistance = fA11+(2.0)*fB1+fC;
                }
                else if (fB1 >= 0.0)
                {
                    fT = 0.0;
                    fSqrDistance = fC;
                }
                else
                {
                    fT = -fB1/fA11;
                    fSqrDistance = fB1*fT+fC;
                }
            }
        }
        else if (fT < 0.0)  // region 6
        {
            fTmp0 = fA01 + fB1;
            fTmp1 = fA00 + fB0;
            if (fTmp1 > fTmp0)
            {
                fNumer = fTmp1 - fTmp0;
                fDenom = fA00-(2.0)*fA01+fA11;
                if (fNumer >= fDenom)
                {
                    fT = 1.0;
                    fS = 0.0;
                    fSqrDistance = fA11+(2.0)*fB1+fC;
                }
                else
                {
                    fT = fNumer/fDenom;
                    fS = 1.0 - fT;
                    fSqrDistance = fS*(fA00*fS+fA01*fT+(2.0)*fB0) +
                        fT*(fA01*fS+fA11*fT+(2.0)*fB1)+fC;
                }
            }
            else
            {
                fT = 0.0;
                if (fTmp1 <= 0.0)
                {
                    fS = 1.0;
                    fSqrDistance = fA00+(2.0)*fB0+fC;
                }
                else if (fB0 >= 0.0)
                {
                    fS = 0.0;
                    fSqrDistance = fC;
                }
                else
                {
                    fS = -fB0/fA00;
                    fSqrDistance = fB0*fS+fC;
                }
            }
        }
        else  // region 1
        {
            fNumer = fA11 + fB1 - fA01 - fB0;
            if (fNumer <= 0.0)
            {
                fS = 0.0;
                fT = 1.0;
                fSqrDistance = fA11+(2.0)*fB1+fC;
            }
            else
            {
                fDenom = fA00-2.0f*fA01+fA11;
                if (fNumer >= fDenom)
                {
                    fS = 1.0;
                    fT = 0.0;
                    fSqrDistance = fA00+(2.0)*fB0+fC;
                }
                else
                {
                    fS = fNumer/fDenom;
                    fT = 1.0 - fS;
                    fSqrDistance = fS*(fA00*fS+fA01*fT+(2.0)*fB0) +
                        fT*(fA01*fS+fA11*fT+(2.0)*fB1)+fC;
                }
            }
        }
    }

    // account for numerical round-off error
    if (fSqrDistance < 0.0)
    {
        fSqrDistance = 0.0;
    }

    //m_kClosestPoint0 = m_rkVector;
    *nearest = v0+ fS*kEdge0 + fT*kEdge1;
    //m_afTriangleBary[1] = fS;
    //m_afTriangleBary[2] = fT;
    //m_afTriangleBary[0] = 1.0 - fS - fT;

    return fSqrDistance;
}

#endif
