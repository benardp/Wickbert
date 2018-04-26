/**
 * @file MLSSurface.h
 *
 * @brief Description of MLSSurface class and supporting structures.
 * @author Matei Stroila based on Andrew Colombi
 */

#ifndef MLSSURFACE_H
#define MLSSURFACE_H

#include <vector>
#include <gsl_matrix.h>
#include <gsl_eigen.h>
#include <gsl_cblas.h>
#include <gsl_blas.h>
#include <gsl_errno.h>
#include <gsl_math.h>
#include <gsl_min.h>
#include <sstream>
#include <fstream>
#include <sfl.h>

#include "Constants.h"
#include "SpaceDef.h"
#include "Kernel.h"
#include "Surface/KDTree.h"

#include <GL/gl.h>

#include "Surface/Implicit/Implicit.h"

#define FUDGE 0.00001       ///< Fudge factor for computing best fit planes
#define SURFACE_TOL 0.00001 ///< "Converged" threshold
#define PROJ_STEPS 10       ///< Number of projection steps before we give up



namespace MLS {

    //class prototype
    class MLSSurface;

    /**
    * @brief Represents a "Surflet" sample.
    *
    * PointSample represents a Surflet sample, i.e. a point with a normal.
    */
    class PointSample : public Point {
        
        public:
            Vector n; //The normal of the sample
        
            PointSample(Point p, Vector m) : Point(p), n(m) {}
            PointSample() : Point(), n() {}
            ~PointSample() {}
    };

    /**
    * @brief Used to interface to GSL.
    *
    * This structure represents the parameters to a function that
    * GSL will seek minima over.
    */
    struct ParamStruct {
        
        MLSSurface *surface;
        double (MLSSurface::*f)(const Point&, Vector&, Kernel*);
        Point q;
        Vector n;
        Kernel *k;
    };

    /**
    * @brief Used to interface to GSL.
    * 
    * This interface represents a function that GSL will
    * seek a minima over.
    */
    double gslMinFuncInterface(double x, void *param);

    /**
    * @brief MLSSurface is responsible for computing an MLS Surface and
    * representing in Implicitly.
    *
    * MLSSurface contains methods for computing the MLS Surface of point set
    * data (typically with normals).  The surface is defined as the set of stationary
    * points of a projection operation.  The projection operation is an
    * optimization process which seeks to project a point onto the surface
    * by repeatedly seeking out a local minimum of an energy function.  The
    * domain of the energy function is restricted to a line defined by a vector
    * field. Several energy and vector field functions are defined in this class.
    * See function documentation for further details.
    */
	class MLSSurface : public Implicit {

    public:
        typedef PointSample Sample;                 ///< Typedef to let use change what kind of samples we use
        typedef KDTree<DIMENSION, Point> KDTree;    ///< KDTree for computing k-nearest queries.
        typedef MLSSurface SurfaceType;

    public:
        
        //construction & destruction
        MLSSurface();
        virtual ~MLSSurface();

        //file reading
        virtual bool readImplicit(std::ifstream &file, bool verbose);
        void readMLSFile(std::ifstream &file);

#ifdef WB_USE_SFL   //only reads SFL if we have it installed
        void readSFLFile(const char *psFile);
        void readSFLSurfletSet(sfl::InStream *in, std::vector<Point> &samplesAsPts, double &scale);
        void scaleSFLInput(std::vector<Point> &samplesAsPts, double scale);
#endif

        //computes the MLS energy
        double computeMLSEnergy(const Point &query,Vector &n,Kernel *kernel);
        double computeMAHEnergy(const Point &query,Vector &n,Kernel *kernel);

        //computes vector fields
        bool computeCovNormal(const Point &query,Kernel *kernel, Vector &normal, Point &avg);
        bool computeOptMLSDir(const Point &query,Kernel *kernel, Vector &normal, Point &avg);
        bool computeAvgNormal(const Point &query,Kernel *kernel, Vector &normal, Point &avg);

        //computes distance from sample points
        bool computeDistanceAvg(const Point &query, Kernel *kernel, double &distAvg);

        //energy minimation processes
        bool projectPoint(const Point &query, Kernel *kernel, bool (MLSSurface::*optDir)(const Point&, Kernel*, Vector& ,Point&),
            double (MLSSurface::*energy)(const Point&, Vector&, Kernel*), Point &proj); 
        double findEnergyMin(const Point &query, Kernel *kernel, Vector &normal,double bound, 
            double (MLSSurface::*energy)(const Point&, Vector&, Kernel*)); 

        //data drawer
        void renderPoints(float *color);

        //gets the kernel
        Kernel *getKernel();

    protected:

		//callback to load file with surfels
		class LoadButtonCallback : public SurfParamButton::Callback
		{
		public:
			MLSSurface *mlsSurf;
			LoadButtonCallback(MLSSurface *pmlsSurf) {mlsSurf = pmlsSurf;}
			virtual void onbuttonpress() {mlsSurf->load();}
		};

		//callback to apply the parameters' changes
		class ApplyParameterChangesCallback : public SurfParamButton::Callback
		{
		public:
			MLSSurface *mlsSurf;
			ApplyParameterChangesCallback(MLSSurface *pmlsSurf) {mlsSurf = pmlsSurf;}
			virtual void onbuttonpress() {mlsSurf->applyParameterChanges();}
		};

		class KernelSelector : public SurfParamComboBox::Callback
		{
		protected:
			Kernel ** _basicKernel;
		public:
			KernelSelector(Surface * surface, Kernel** rf);
			virtual ~KernelSelector();
			virtual void itemselected(unsigned int selection);
		};


		//functions for callbacks:
		void load();
		void applyParameterChanges();

        int m_VertexNr;         ///< Number of input vertices

        KDTree *m_Tree;         ///< The kd-tree for searching the point samples
        
        std::vector<Sample> m_Samples;  ///< All the point samples in a vector
        
        bool m_InputNormals;    ///< Flag to indicate whether we have input normals

        Kernel *m_Kernel;        ///< Support kernel for the MLS Surface.

		std::string filename; ///< SFL input file name

		bool kernelIsFinite; ///< flag finite or not kernel
		int kernelMaxNumberOfNeighbors; ///< max # of neighbors
		double kernelSqrSize;  ///< the squared size of a finite kernel
		double kernelSqrVariance; ///< the squared variance of a Gauss kernel
		double kernelLinearCutOff; ///< When the linear clipper kicks in (asbolute distance from 0)
    };
};



#endif // MLSSURFACE_H
