/**
 * Implementation of a volume data implicit 
 * @file ITKImplicit.h
 * @date 10/01/2005
 * @author Matei N. Stroila modified a lot by Elmar Eisemann :)
 */

#ifndef ITKImplicit_h
#define ITKImplicit_h

#include "Surface/Implicit/Implicit.h"

#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkImageSpatialObject.h" 

typedef double PixelType;
const unsigned int      Dimension = 3;
typedef itk::Image< PixelType, Dimension >      ImageType;
typedef itk::ImageSeriesReader< ImageType >     ReaderType;
typedef itk::BSplineInterpolateImageFunction<ImageType, double> InterpFunc;
typedef itk::GDCMImageIO                        ImageIOType;
typedef itk::GDCMSeriesFileNames                NamesGeneratorType;
typedef itk::ImageSpatialObject<Dimension ,PixelType> ImageSpatialObjectType; 

class ParticleBoundingBox;
class Implicit;
class SurfaceMesh;

class ITKImplicit : public Implicit
{
public:
	ITKImplicit();
	
	double proc(const gmVector3 & x);
	gmVector3 grad(const gmVector3 & x);
	gmMatrix3 hess(const gmVector3 & x);

	//Initialize the implicit
	void init(void);
	void initFromImplicit(Implicit * imp);
	void initFromSurfaceMesh(SurfaceMesh * mesh);

	void applyParameterChanges();
	
	MAKE_NAME();
	
private:
	//callback to initialize the implicit
	class SetITKImplicit : public SurfParamButton::Callback
	{
	public:
		ITKImplicit *itkImp;
		SetITKImplicit(ITKImplicit *pitkImp) {itkImp = pitkImp;}
		virtual void onbuttonpress() {itkImp->init();}
	};

	//callback to apply the affine transformation the implicit
	class ApplyParameterChangesCallback : public SurfParamButton::Callback
	{
	public:
		ITKImplicit *itkImp;
		ApplyParameterChangesCallback(ITKImplicit *pitkImp) {itkImp = pitkImp;}
		virtual void onbuttonpress() {itkImp->applyParameterChanges();}
	};

	//callback to apply the affine transformation the implicit
	class LoadButtonCallback : public SurfParamButton::Callback
	{
	public:
		ITKImplicit *itkImp;
		LoadButtonCallback(ITKImplicit *pitkImp) {itkImp = pitkImp;}
		virtual void onbuttonpress() {itkImp->load();}
	};

		//callback to apply the affine transformation the implicit
	class SaveButtonCallback : public SurfParamButton::Callback
	{
	public:
		ITKImplicit *itkImp;
		SaveButtonCallback(ITKImplicit *pitkImp) {itkImp = pitkImp;}
		virtual void onbuttonpress() {itkImp->save();}
	};


	//initialize myImage such that it fits the bounding box (min, max)
	void initScalingAndOriginFromBoundingBox(const gmVector3 & min, const gmVector3 & max);

	//load a file
	//save a file
	void load();
	bool loadDicom(const std::string & dicomDir);
	void save() const;

	
	//Dicom directory
	std::string dirFile;
	//ITK reader
	ReaderType::Pointer reader;
	//DICOM image factory
	ImageIOType::Pointer dicomIO; 
	//pointer to the volume data
	ImageType::Pointer myImage;
	//B-spline interpolation function
	InterpFunc::Pointer myInterpFunc;
	//Spatial object for the volume data
	ImageSpatialObjectType::Pointer imageSO;
	//Bounding box of the volume data
	ImageSpatialObjectType::BoundingBoxType *boundingBox; 
	//Used for reading the DICOM data
	NamesGeneratorType::Pointer nameGenerator;
	
	//CT number
	double threshold;
	//Volume scale factor
	double scale;
	//Volume translation
	gmVector3 tr;
	//Spline order (default 3)
	int order;
	//Particle attribute to set the bounding box of the particle system 
	//equal to the bounding box of the volume data.
	ParticleBoundingBox *bbox;

	//SurfaceReference, in the case that no filename is specified
	//the ITK implicit is generated from the implicit
	Surface * surface;

	//modified center after affine transformation
	//storing this value speeds up a little.
	gmVector3 center;

	ImageType::PointType originalOrigin;
	ImageType::SpacingType originalSpacing;

	//wanted resoltion of the voxel grid
	double boundingBoxScaling;
	int voxelRes[3];
};

#endif 

