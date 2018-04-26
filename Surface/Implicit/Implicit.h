/**
 * @file Implicit.h
 * @author Nate Carr, Jeff Decker, John C. Hart, and Terry Fleury
 * Implicit Class definition
 */

/** @page implicitpage The Implicit Object Library
 *
 * The implicit class is a class that represents functions.
 * This class is typically invoked as "Implicit *f" where f is a pointer to
 * some instance of a concrete subclass of Implicit. The function can be
 * evaluated as f->proc(x) or *f(x) where x is a point in space.  In this
 * implementation, we define the value of the function to be negative
 * "within" an implicit surface, 0 "on" the surface,  and positive
 * "outside" the surface.
 *
 * Derivatives grad(x) and hess(x) are also available. Some implicits define
 * these analytically, or they may otherwise default to a slower numerical
 * version.
 *
 * The parameters of each Implicit object are available through a uniform
 * parameter interface:
 * @li int qlen() - the number of parameters
 * @li setq(double *q) - set the parameters from the array q.  This also
 *     saves a copy of the previous q array for evaluation of the curve over
 *     "time".  External classes will call this method.
 * @li _setq(double *q) - set the parameters from the array q.  This version
 *     does not save a copy of the previous q array and should be overridden
 *     by most subclasses.
 * @li getq(double *q) - return current parameters in the array q. You are
 *     responsible for allocating q. Make sure it is at least qlen() long.
 * @li getqold(double *q) - return the previous parameters of the surface in
 *     the array q.  You must allocate space for the q array.
 * @li getqname(char **qn) - load an array of strings with names of the
 *     parameters. You are responsible for allocating an array of (char *)
 *     in qn[]. The Implicit object will allocate the individual strings.
 * @li procq(gmVectro3 x, double *dq) - returns in dq an array of length
 *     qlen() the gradient of proc() wrt its parameters, evaluated at the
 *     point x.
 *
 * Interval extensions proc(Box<double> b), grad(Box<double> b) and 
 * hess(Box<double> b) (and eventually procq()) should also be available for
 * Implicit objects.  Additionally, if you are planning on doing 4D critical
 * point finding, you need to implement procqt() = dproc(X,t)/dq(t) and 
 * gradqt =  dgrad(X,t)/dq(t).  
 *
 * Some implicit objects are operators on other implicits. The member
 * numChildren() returns the number of children the operator currently has
 * defined, and maxChildren() returns the number of children the operator
 * can handle (-1 if it is unlimited). Children are accessed with
 * getChild()/setChild()/getChildren()/setChildren().
 *
 * Objects also have a name associated with them, and the name is
 * synthesized uniquely. The name is available via getObjectName() and its
 * variations.
 *
 * Information about creating a new Implicit object is available in the
 * detailed documentation of the Implicit class.
 *
 * @note As of February 3, 2003, you have the option of using a new calling
 * hierarchy for proc, grad, and hess.  Previously, you had to define proc,
 * grad, and hess for both gmVector3 and Box3d (and Box4d).  If you #define
 * INTERVAL_EVAL_ONLY, all you need to define is proc, grad, and hess for
 * Box<double> and the other (non-interval) methods use these.  For example,
 * proc(gmVector3) calls proc(Box<double>) by converting the gmVector to a
 * Box<double>, calling proc(Box<double>), and then returning the center()
 * of the resulting Intervald.  If you undefine the compiler option
 * INTERVAL_EVAL_ONLY, you use the old calling hierarchy (assuming that all
 * of the gmVector3 versions of proc, grad, and hess have been maintained to
 * keep equivalent with the Box<double> versions).
 */


#ifndef __IMPLICIT_H
#define __IMPLICIT_H

#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include <string>
#include <valarray>
#include <vector>
#include <sstream>
#include <fstream>

#include "tools/factory.h"
#include "pixmaps.h"
#include "Surface/Surface.h"

/**
 * If you are using the standard versions of grad() and hess() (ie. you have
 * not defined INTERVAL_EVAL_ONLY), you may #undef VIRTUAL_DERIV to force
 * grad() and hess() to use numerical approximations of the functions EVEN IF
 * these methods are overridden in subclasses.  You may want to do this
 * when testing out new implemenations of grad() and hess() to verify that
 * they are accurate.
 */
#define VIRTUAL_DERIV virtual /* Allow sublass to override grad and hess */
// #undef VIRTUAL_DERIV       /* Force numerical approx of grad and hess */

typedef std::valarray<double>       DoubleArray;
typedef std::vector<double>         DoubleVector;
typedef std::vector< std::string >  NameVector;

/** The master class of all implicit surfaces.
 *
 * The implicit class is a class that represents functions.
 * This class is typically invoked as "Implicit *f" where f is a pointer to
 * some instance of a concrete subclass of Implicit. The function can be
 * evaluated as f->proc(x) or *f(x) where x is a point in space.  In this
 * implementation, we define the value of the function to be negative
 * "within" an implicit surface, 0 "on" the surface,  and positive
 * "outside" the surface.
 *
 * Subclasses of Implicit need to define the function defining the implicit
 * surface: double proc(const gmVector3 & x) and should define an interval
 * proc(Box<double> b) as well.
 * 
 * Subclasses also need to maintain their own parameterization, and so need
 * to define the following parameterization interface:
 * @li unsigned int qlen() - the number of parameters
 * @li setq(double *) - set the parameters from the given array of double
 * @li getq(double *) - load an array of double with current parameters
 * @li getqname(char **) - load an array of strings with names of the parameters
 * @li procq(gmVectro3 x, double *dq) - returns in dq an array of length
 *     qlen() the gradient of proc() wrt its parameters, evaluated at the
 *     point x.
 * 
 * Subclasses should define their own grad(), hess() and procq() methods,
 * though numerical approximations are implemented in Implicit.
 *
 * \nosubgrouping
 */
class Implicit : public Surface 
{   
  public:
    /// @name Constructors and Destructor
    /// @{
    Implicit();     
    virtual ~Implicit();
    /// @}

	//added this is for passing events from a GUI -ms
	virtual void processEvent(const std::string&){};
	
	/** Reads implicit parameters from a file.
	 * The default implementation grabs first the number of parameters
	 * and then the individual parameters (as doubles). Child classes
	 * may override this to provide any kind of detailed loading capabilities
	 * necessary. This is called immediately after the class's name is
	 * discovered in the file, and the class is instantiated. Implementations
	 * of this function should only consume exactly the amount of input
	 * necessary. Implicit hierarchy construction is still handled in
	 * ImpFileManager::readImplicit()
	 * \param file The input file stream.
	 * \param verbose Whether or not the function should output error messages to stderr.
	 * \return Whether or not the load succeeded.
	 * \author Scott Kircher
	 */
	virtual bool readImplicit(std::ifstream &file,bool verbose);

	/** Common param loading.
	 * Many child classes may need to have initial behavior similar to Implicit::readImplicit()
	 * and can use this utility function to read an array of doubles from the file
	 * they may then do whatever they wish with it (readImplicit() just calls setq using this array).
	 * This function reads first the number of parameters to load, and then the actual paramters.
	 * \param file The input file stream.
	 * \param params Vector of doubles that will be filled with the paramters read from the file.
	 * \author Scott Kircher
	 */
	void readParameters(std::ifstream &file,std::vector<double> &params);

    /** @name Evaluation
     * An implicit is evaluated by calling its proc(x) method which returns
     * the scalar value associated with 3-D point x. The implicit object
     * itself can be called as a function, via the overloaded parentheses
     * operator.
     * 
     * The function proc() can also be differentiated with respect to x
     * once: grad() or twice: hess(), or with respect to its parameters:
     * procq(). The functions proc() and its derivative grad() can also be
     * differentiated with respect to the Qs over time as procqt() and
     * gradqt() (assuming a linear trajectory between the last two parameter
     * settings over unit time).
     *  @{ */

	virtual double proc(const gmVector3 & x) = 0;
    virtual Intervald proc(const Box<double>& b);

    double operator()(const gmVector3 & x);
    VIRTUAL_DERIV gmVector3 grad(const gmVector3 & x);
    VIRTUAL_DERIV gmMatrix3 hess(const gmVector3 & x);

    Intervald operator()(const Box<double>&  b);
    virtual Box3d grad(const Box<double>& x);
    virtual IMatrix3d hess(const Box<double>& x);
    //I don't think this function is general enough to stay here.
	//It calculates the hessian, and the grad (even using Fx Fy Fz)
	//it should be avoided to add all kind of functions like this to 
	//the implict. It could nevertheless be added to ImplicitInterrogator
	//but even there it lacks some generality... - Elmar
	//gmVector3 gradGradN2(const gmVector3 &x);
    /// @}

    /// @name Differentiation @{
    virtual void procq(const gmVector3&, double*); ///< Deriv of func wrt parameters
    void procq(const gmVector3&, DoubleArray&);    ///< Evaluation of dFdQ into array
    void procq(const gmVector3&, DoubleVector&);   ///< Evaluation of dFdQ into vector
    void procq(const gmVector3&, TNT::Vector<double>&);  ///< Eval of dFdQ into TNT

    virtual void procq(const Box<double>& , Intervald*); 
    virtual void gradq(const Box<double>& , Intervald*, Intervald*, Intervald*);

    /** Derivative of function wrt Qs over time.
     * Assumes a linear interpolation between parameters qold and q over
     * one unit of time.  No need to override this method.  Instead,
     * override procq().
     */
    virtual Intervald proct(const Box<double>& );

    /** Derivative of gradient wrt Qs over time.
     * Assumes a linear interpolation between parameters qold and q
     * over one unit of time.  No need to override this method.  Instead,
     * override gradq().
     */
    virtual Box3d gradt(const Box<double>& );

    /*! Returns the value of epsilon used for numerical derivative approximation.
     *  \return m_epsilon
     *  \author Jared Hoberock
     */
    inline double getEpsilon(void) const
    {
      return m_epsilon;
    } // end Implicit::getEpsilon()

    /*! This method takes the partial deriviative of Hess wrt variable i.
     *  \param x A point in space.
     *  \param i A variable to differentiate wrt. x: 0, y: 1, z: 2.  For
     *           example, Hess_x f(x) would be hessi(x,0), etc.
     *  \return The partial derivative of the hessian at x wrt variable i.
     *  \author Jared Hoberock
     */
    gmMatrix3 hessi(const gmVector3 &x, const unsigned int i);

  protected:
    virtual void dqtdt(double *);
    /// @}

    /** @name Differential Geometry
     * These members provide a variety of numerical and analytic
     * differential geometry operations and parameters. They are all defined
     * via proc() and its derivatives.
     * @{ */
  public:
    virtual double lipschitz(const Box<double>& x);
    void setEpsilon(double epsilon);
    gmVector3 normal(const gmVector3 & x);
    gmVector3 tangent(const gmVector3 & x);
    gmVector3 binormal(const gmVector3 & x);
    double k(const gmVector3 & x);
    double gaussianCurvature(const gmVector3 & x);
	Intervald gaussianCurvature(const Box<double>& x);
	//for gaussian curvature zero sets
	double numeratorGaussianCurvature(const gmVector3 & x);
	Intervald numeratorGaussianCurvature(const Box<double>& x);

    /*! This method approximates the gradient of gaussian curvature.
     *  \param x A point in space.
     *  \return An approximation of the gradient of gaussian curvature at x.
     *  \author Jared Hoberock
     */
    gmVector3 gradGaussianCurvature(const gmVector3 & x);

    double meanCurvature(const gmVector3 & x);

    /*! This method approximates the gradient of mean curvature.
     *  \param x A point in space.
     *  \return An approximation of the gradient of mean curvature at x.
     *  \author Jared Hoberock
     */
    gmVector3 gradMeanCurvature(const gmVector3 & x);

    double normalCurvature(const gmVector3 & x, const gmVector3 & v);

    /*! This method approximates the gradient of normal curvature.
     *  \param x A point in space.
     *  \param v A vector in the tangent space of this Implicit.
     *  \return An approximation of the gradient of normal curvature at x in direction v.
     *  \author Jared Hoberock
     */
    gmVector3 gradNormalCurvature(const gmVector3 & x, const gmVector3 & v);

    gmVector3 shapeOperator(const gmVector3 & x, const gmVector3 & v);

	 /*! This method computes the radial curvature.
     *  \param x A point in space.
     *  \param v View vector.
     *  \return The radial curvature at x in direction v.
     *  \author Matei Stroila
     */
	double radialCurvature(const gmVector3 & x, const gmVector3 & v);
	double numeratorOfRadialCurvature(const gmVector3 & x, const gmVector3 & v);
	 /*! This method computes the gradient of the radial curvature.
	 *  \param gradRadCurv The gradient of the radial curvature at x in direction v.
     *  \param x A point in space.
     *  \param v View vector.
     *  \author Matei Stroila
     */
    void gradRadialCurvature(gmVector3& gradRadCurv, const gmVector3 & x, const gmVector3 & v);
	void gradNumeratorOfRadialCurvature(gmVector3& gradRadCurv, const gmVector3 & x, const gmVector3 & v);

    /// @}

    /** @name Parameters
     * Provides a uniform interface to the parameters of all Implicit objects.
     * Parameters are always maintained as an array of floating point values.
     * @{ */

    /** Number of implicit model parameters for this implicit object including
     * its children.
     */
    virtual unsigned int qlen();

    /** Number of implicit model parameters for this implicit object only,
     * not including its children.
     */
    virtual unsigned int plen();

    virtual void getq(double*);  ///< Get copy of the implicit model parameters
    void getq(DoubleArray&);     ///< Get parameters into a double array
    void getq(DoubleVector&);    ///< Get parameters into a double vector
    void getq(TNT::Vector<double>&); ///< Get parameters into TNT double vector

    /** The previous implicit parameters. The member variable qold gets
     *  updated everytime setq is called with the existing parameters that
     *  have been replaced by the new parameters from the setq call.
     */
  private:
    double* qold;                  ///< Array to hold "previous" Qs
    bool qsChangedFlag;            ///< Set to true when setq is called
    int qoldsize;                  ///< Size of qold "array"
    void allocateQold(void);       ///< Get memory space for qold "array"
  public:
    virtual void getqold(double*); ///< Get the OLD implicit model parameters
    void getqold(DoubleArray&);    ///< Get OLD parameters into a double array
    void getqold(DoubleVector&);   ///< Get OLD parameters into a double vector
    void getqold(TNT::Vector<double>&); ///< Get OLD parameters into TNT vector
    bool qsChanged(void);          ///< Is qold different from current Qs?

    /** The difference between _setq() and setq() is that setq() SAVES the
     *  previous set of Qs into the qold array, but _setq() does NOT.
     *  _setq() is called directly only in subclasses of Implicit.  Outside
     *  the Implicit class hierarchy, always call setq().  So, if you
     *  want to save the old set of Qs, a subclass should override _setq().
     *  In fact, the default _setq() in Implicit does nothing, so you have
     *  to override it to do the work of setting the current Q array.
     *  However, on the off chance you have a subclass where you do not want
     *  to save the old Qs, then override setq() to do the work of settings
     *  the current Qs.  Typically you will want to save the old Qs so
     *  subclasses will usually override (and call) _setq().
     */
    virtual void _setq(double*); ///< Overridden by subclasses to set Qs
    virtual void setq(double*);  ///< Change the implicit model parameter vector
    void setq(DoubleArray&);     ///< Set parameters from a double array
    void setq(DoubleVector&);    ///< Set parameters from a double vector
    void setq(TNT::Vector<double>&); ///< Set parameters from TNT double vector

    virtual void getqname(char** qn); /// The name of each parameter
    void getqname(NameVector&);       ///< Fetches parameter names into a vector

	/** Change the implicit surface according to the control points of its
	 * formulation. This defaults to a version that modifies the q-params of
	 * an implicit to pass through "controller" particles. However, for RBF's
	 * this method solves the linear system to compute the RBF weights.
	 *
	 * floaters query and display implicit through ImplicitInterrogators
	 * controllers call interpolate via SurfaceAdhesion through ImplicitInterrogator
	 *
	 * centers call interpolate via SurfaceAdhesion through ImplicitInterrogator
	 */
	virtual void interpolate(Particles *ps, float phi);

    /// @}

    /** @name Hierarchy
     * Implicit objects are often composites of other Implicit objects.
     * For example, an Implicit might be the blended union of other Implicit
     * primitives.
     * @{ */
    virtual int numChildren(); ///< Number of children currently defined
    virtual int maxChildren(); ///< Max number of children that can be defined
    virtual bool setChildren(std::vector<Implicit*>); ///< If Implicit has any
    /// If Implicit has children, sets the child specified by index.
    virtual bool setChild(int index, Implicit* child);
    virtual void getChildren(std::vector<Implicit*>&);///< If Implicit has any
    virtual Implicit* getChild(int index);  ///< Get the specified child

    /// Returns the Pixmap for the Surface as a char* in the XPM image format.
    virtual const char** getPixmapXPM(const int& size) const;

	/** Holds the color of the implicit object. If alpha is less than one
     * then the color should be alpha*RGB + (1-alpha)*(combination of RGB of
     * children).
     */
    gmVector4 color;
    /// @}

    ///// @name Interrogation
    ///// @{
    ///// Polygonizes the Implicit into the mesh stored in Surface::m_mesh
    //char* polygonize(double size, int bounds, gmVector3 x, 
    //                 int (*triproc)(int,int,int,jbVertices));
 
    ///// Polygonizes specifically using marching cubes
    //char* marchingcubes(double size, gmVector3 x0, gmVector3 x1, 
    //                    int (*triproc)(int,int,int,jbVertices));

    /** 
     * Determines the surface area of the Implicit. For now, -1 means area
     * cannot be calculated. 
     * @returns -1 if area cannot be calculated
     * @todo implement this
     */
    double area() { return -1.0; }
    /// @}
  public:

    /// @name Derivative Helpers
    /// @{

    /** Used for numerical derivative approximations.
     * Also used for subclass derivative numerical approxations e.g. blend.
     *  @see grad, hess
     */
    double m_epsilon;

    double Fx(const gmVector3 & x);
    double Fy(const gmVector3 & x);
    double Fz(const gmVector3 & x);

    double Fxx(const gmVector3 & x);
    double Fxy(const gmVector3 & x);
    double Fxz(const gmVector3 & x);

    double Fyx(const gmVector3 & x);
    double Fyy(const gmVector3 & x);
    double Fyz(const gmVector3 & x);

    double Fzx(const gmVector3 & x);
    double Fzy(const gmVector3 & x);
    double Fzz(const gmVector3 & x);

    /*! This method takes the partial derivative of Fxx wrt variable i.
     *  \param x A point in space.
     *  \param i A variable to differentiate wrt. x: 0, y: 1, z: 2.
     *  \return The partial derivative of Fxx wrt variable i.
     */
    double Fxxi(const gmVector3 & x, unsigned int i);

    /*! This method takes the partial derivative of Fxy wrt variable i.
     *  \param x A point in space.
     *  \param i A variable to differentiate wrt. x: 0, y: 1, z: 2.
     *  \return The partial derivative of Fxy wrt variable i.
     */
    double Fxyi(const gmVector3 & x, unsigned int i);

    /*! This method takes the partial derivative of Fxz wrt variable i.
     *  \param x A point in space.
     *  \param i A variable to differentiate wrt. x: 0, y: 1, z: 2.
     *  \return The partial derivative of Fxz wrt variable i.
     */
    double Fxzi(const gmVector3 & x, unsigned int i);

    /*! This method takes the partial derivative of Fyx wrt variable i.
     *  \param x A point in space.
     *  \param i A variable to differentiate wrt. x: 0, y: 1, z: 2.
     *  \return The partial derivative of Fyx wrt variable i.
     */
    double Fyxi(const gmVector3 & x, unsigned int i);

    /*! This method takes the partial derivative of Fyy wrt variable i.
     *  \param x A point in space.
     *  \param i A variable to differentiate wrt. x: 0, y: 1, z: 2.
     *  \return The partial derivative of Fyy wrt variable i.
     */
    double Fyyi(const gmVector3 & x, unsigned int i);

    /*! This method takes the partial derivative of Fyz wrt variable i.
     *  \param x A point in space.
     *  \param i A variable to differentiate wrt. x: 0, y: 1, z: 2.
     *  \return The partial derivative of Fxy wrt variable i.
     */
    double Fyzi(const gmVector3 & x, unsigned int i);

    /*! This method takes the partial derivative of Fzx wrt variable i.
     *  \param x A point in space.
     *  \param i A variable to differentiate wrt. x: 0, y: 1, z: 2.
     *  \return The partial derivative of Fzx wrt variable i.
     */
    double Fzxi(const gmVector3 & x, unsigned int i);

    /*! This method takes the partial derivative of Fzy wrt variable i.
     *  \param x A point in space.
     *  \param i A variable to differentiate wrt. x: 0, y: 1, z: 2.
     *  \return The partial derivative of Fzy wrt variable i.
     */
    double Fzyi(const gmVector3 & x, unsigned int i);

    /*! This method takes the partial derivative of Fzz wrt variable i.
     *  \param x A point in space.
     *  \param i A variable to differentiate wrt. x: 0, y: 1, z: 2.
     *  \return The partial derivative of Fzz wrt variable i.
     */
    double Fzzi(const gmVector3 & x, unsigned int i);

    Intervald Fx(const Box<double>& );
    Intervald Fy(const Box<double>& );
    Intervald Fz(const Box<double>& );

    Intervald Fxx(const Box<double>& );
    Intervald Fxy(const Box<double>& );
    Intervald Fxz(const Box<double>& );

    Intervald Fyx(const Box<double>& );
    Intervald Fyy(const Box<double>& );
    Intervald Fyz(const Box<double>& );

    Intervald Fzx(const Box<double>& );
    Intervald Fzy(const Box<double>& );
    Intervald Fzz(const Box<double>& );
	
	    //@}

public:
		
	Intervald Fxxi(const Box<double>& , unsigned int i);
    Intervald Fxyi(const Box<double>& , unsigned int i);
    Intervald Fxzi(const Box<double>& , unsigned int i);
	
    Intervald Fyxi(const Box<double>& , unsigned int i);
    Intervald Fyyi(const Box<double>& , unsigned int i);
    Intervald Fyzi(const Box<double>& , unsigned int i);
	
    Intervald Fzxi(const Box<double>& , unsigned int i);
    Intervald Fzyi(const Box<double>& , unsigned int i);
    Intervald Fzzi(const Box<double>& , unsigned int i);
	
	Intervald D3F(const Box<double>& x, unsigned int i, unsigned int j, unsigned int k);
	double D3F( const gmVector3& x,  unsigned int i,  unsigned int j,  unsigned int k) ;
	
	virtual void output(std::ostream &out);
	virtual void input(std::istream &in);
};

/** Computes the outer product between to 3-vectors. The result might
 * be the transpose of what it really should be. I haven't checked the
 * ordering.
 * @param x 3-vector factor
 * @param y 3-vector factor
 * @return 3x3 matrix containing pairwise products of x,y elements
 * @note Hate to put this here of all places, but want to keep libgm unchanged
 * and this is the only file I know will be included by everyone.
 * ELMAR: err... not by me... why don't we at least move it up as much as possible in the hierarchy?
 * Surface.h, is also the first place where the vectors from libgm are included...
 * I put it there with the comment from above... once it is accepted by everyone we can delete this
 * actually if libgm is extended, we could just do this in a separate libgmext.h
 * which could then be included instead of libgm.
 */
//gmMatrix3 outer(const gmVector3 & x, gmVector3 y);
/// Interval version of outer product of two Boxes
//IMatrix3d outer(Box3d x, Box3d y);

/** Shorthand type for a factory that creates Implicits.
 * ImplicitFactory is a genericFactory<Implicit>.
 */
typedef genericFactory<Implicit> ImplicitFactory;

/** Shorthand type for the ImplicitFactory registry.
 * This type is essentially a map of classname strings to
 * functions that create the new class.
 * \sa IMPLICIT_REGISTRY
 */
typedef genericFactory<Implicit>::FN_REGISTRY ImplicitRegistry;

/** Macro to register a new Implicit subclass.
 * \param class_name The token name of the class.
 * \param class_string The string name of the class.
 * \note This belongs in the .cpp file, right after the includes.
 */
#define REGISTER_IMPLICIT(class_name,class_string) \
  int grab##class_name; \
    std::string class_name::registry_name = class_string; \
  namespace { \
  registerInFactory<Implicit,class_name> registerMe(class_name::registry_name);\
  }

/** Macro to add name() and objectName functionality to all subclasses of
 * Implicit.
 */
#define MAKE_NAME() \
  static std::string registry_name; \
  virtual const std::string name() { return registry_name; } \
  virtual void resetObjectName() { \
  std::ostringstream objectNumberName; \
  objectNumberName << ":" << ++objectNumber; \
  objectName = Implicit::defaultObjectNamePrefix + name() + \
                     objectNumberName.str(); \
  }

/** Macro to get the registry for the ImplicitFactory.
 * The registry is a map of classname strings to new class constructors.
 * \return The registry of type ImplicitRegistry.
 */
#define IMPLICIT_REGISTRY genericFactory<Implicit>::instance().registry

/** Macro to construct a new subclass of Implicit.
 * \param class_string String name of the class to create.
 * \return A std::auto_ptr<Implicit> to the appropriate subclass.
 * This is pretty much just a Implicit* but the auto_ptr adds some
 * safety mechanisms to avoid memory leaks.
 */
#define NEW_IMPLICIT(class_string) \
  genericFactory<Implicit>::instance().create(class_string);

#endif 

