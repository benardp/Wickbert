/**
 * @file MLSSurface.cpp
 *
 * @brief Implementation of MLSSurface class.
 * @author Matei Stroila based on Andrew Colombi
 */
#ifdef WB_USE_SFL

#include "MLSSurface.h"

using namespace MLS;

/**
 * @brief Creates an MLSSurface instance.
 */
MLSSurface::MLSSurface() {

    m_VertexNr = 0;

    m_InputNormals = false;

    m_Tree = new KDTree();

	new SurfParamString(this,&filename,"","filename", "Input SFL file");
	new SurfParamButton(this,new LoadButtonCallback(this),"load","","load file (extension .sfl suggested)");
	
	new	SurfParamComboBox(this,new KernelSelector(this, &m_Kernel), "KernelSelector","KernelSelector",
		"The selection decides which kernel to use in MLS");

	new SurfParamBool(this,&kernelIsFinite,false,"kernelIsFinite", "Finite Kernel");
	new SurfParamInt(this,&kernelMaxNumberOfNeighbors,15,"kernelMaxNumberOfNeighbors", "Max number of neighbors");
	new SurfParamDouble(this,&kernelSqrSize,0.1,"kernelSqrSize", "Kernel Square Size");
	new SurfParamDouble(this,&kernelSqrVariance,0.25,"kernelSqrVariance", "Gauss Kernel Squared Variance");
	new SurfParamDouble(this,&kernelLinearCutOff,0.1,"kernelLinearCutOff", "GaussClip Kernel linear cutoff");

	new SurfParamButton(this,new ApplyParameterChangesCallback(this),"applyParameterChanges","","Apply Parameters' Changes");


}

/**
 * @brief Frees an MLSSurface instance.
 */
MLSSurface::~MLSSurface()
{ 
	delete m_Tree;
	delete m_Kernel;
}

void MLSSurface::load()
{
	
	readSFLFile(filename.c_str());
}

void MLSSurface::applyParameterChanges()
{
	m_Kernel->setMaxNeighbors(kernelMaxNumberOfNeighbors);
	m_Kernel->setIsFinite(kernelIsFinite);
    m_Kernel->setKernelSize(kernelSqrSize);

	if (GaussKernel *_kernel = dynamic_cast<GaussKernel *>(m_Kernel))
	_kernel->setSqrVariance(kernelSqrVariance);

	if (GaussianClip *_kernel = dynamic_cast<GaussianClip *>(m_Kernel)){
		_kernel->setSqrVariance(kernelSqrVariance);
		_kernel->setLinearCutOff(kernelLinearCutOff);
	}
    
}
/**
 * @brief Reads point data from a file.
 *
 * @param file The file to read the data from.
 * @param verbose Verboose ouput?  Doesn't do anything here.
 *
 * Reads a file containing data points either data points or the name
 * of a file containing a PointShop formatted file. The initial integer (0 or 1) in
 * the file indicates which data it contains.
 */
bool MLSSurface::readImplicit(std::ifstream &file, bool verbose) {

    int mode;
    char psFile[255];

    //The mode of the file tells us whether we read our old format or
    //PointShop files
    file >> mode;
    
    if (mode == 1) {
        
        //HACK!! For the moment this file location is relative to
        //the working directory of the process. !!KCAH
        file >> psFile;
#ifdef WB_USE_SFL
        readSFLFile(psFile);
#else
        std::cerr << "SFL file format not supported in this build." << std::endl;
#endif
    }
    else {
        
        readMLSFile(file);
    }
    
    return true;
}

/**
 * @brief Reads our point data format from a file.
 *
 * @param file The file to read the data from.
 */
void MLSSurface::readMLSFile(std::ifstream &file) {

    int mode;
    int i,j;
    float p;
    Sample curSample;
    std::vector<Point> samplesAsPts;

    //The mode of the file tells us if we read normals
    file >> mode;

    if (mode ==0)
        m_InputNormals = false;
    else
        m_InputNormals = true;
        
    file >> m_VertexNr; //number of points
    
    //clear out any old data
    m_Samples.clear();
    for (i = 0; i < m_VertexNr; i++) {

        for (j = 0; j < DIMENSION; j++) {
            file >> p;
            curSample[j] = p;
        }
        if (m_InputNormals) {
            for (j = 0; j < DIMENSION; j++) {
                file >> p;
                curSample.n[j] = p;
            }
        }
        
        m_Samples.push_back(curSample);
        samplesAsPts.push_back(Point(curSample));
    }

    //clear out old KDTree and make a new one
    if (m_Tree)
        delete m_Tree;
    m_Tree = new KDTree();
    m_Tree->generate(samplesAsPts);
}

#ifdef WB_USE_SFL   //only reads SFL if we have it installed
/**
 * @brief Reads a PointShop formatted file.
 *
 * @param psFile Name of the file containing the data.
 */
void MLSSurface::readSFLFile(const char *psFile) {

    sfl::InStream *in;  // the surfel input stream
    std::vector<Point> samplesAsPts;    // temp storage used to create the KDTree
    double scale;       // data scaling factor to fit into wickberts viewer

    //clear out any old data
    m_Samples.clear();

    if (!(in = sfl::InStream::open(psFile)))
        return;

    in->setReadAsWorldCoordinates(true);

    //read all surfel sets
    int numSurfelSets = in->getNumSurfelSets();
    for (int i = 0; i < numSurfelSets; ++i) {
        
        in->seekSurfelSet(i);
        readSFLSurfletSet(in, samplesAsPts, scale);
    }
    //close the surfel stream
    sfl::InStream::close(in);
    m_VertexNr = m_Samples.size();

    //scale all the input to fit the wickbert bounding volume
    scaleSFLInput(samplesAsPts, scale);

    //clear out old KDTree and make a new one
    if (m_Tree)
        delete m_Tree;
    m_Tree = new KDTree();
    m_Tree->generate(samplesAsPts);
}

/**
 * @brief Reads a single SFL (PointShop) surflet set.
 *
 * @param in Input stream of surflet data.
 * @param samplesAsPts Returned, a vector of the surflet positions.
 * @param scale Scaling factor for all point positions.
 */
void MLSSurface::readSFLSurfletSet(sfl::InStream *in, std::vector<Point> &samplesAsPts, double &scale) {

    Sample curSample;   //temp variable for storing read data.

    in->setSurfelSetPropertyHints(SFLPROP_POSITION | SFLPROP_NORMAL);
    int  rs = in->getSurfelSetRes0Index();

    if (!(in->seekResolution(rs)))
        return;

    int numSurf = in->getResolutionNumSurfels();

    for (int i = 0; i < numSurf; ++i) {

        in->beginSurfel();
        in->readSurfelPosition3(curSample[0], curSample[1], curSample[2]);
        in->readSurfelNormal3(curSample.n[0], curSample.n[1], curSample.n[2]);
        in->endSurfel();

        //update our scaling factor, essential ensure it is the largest coordinate we've
        //seen so far
        scale = (scale > abs(curSample[0])) ? scale : abs(curSample[0]);
        scale = (scale > abs(curSample[1])) ? scale : abs(curSample[1]);
        scale = (scale > abs(curSample[2])) ? scale : abs(curSample[2]);

        m_Samples.push_back(curSample);
        samplesAsPts.push_back(Point(curSample));
    }
}

/**
 * @brief Surfels get scaled by to fit a bounding box of of [1,-1] from [maxP,minP].
 *
 * @param samplesAsPts Samples as points - they get scaled aswell.
 * @param scale Scaling value to apply to all points.
 */
void MLSSurface::scaleSFLInput(std::vector<Point> &samplesAsPts, double scale) {

    scale = 1 / scale;

    for (int i = 0; i < m_VertexNr; ++i) {

        samplesAsPts[i][0] *= scale;
        samplesAsPts[i][1] *= scale;
        samplesAsPts[i][2] *= scale;

        m_Samples[i][0] *= scale;
        m_Samples[i][1] *= scale;
        m_Samples[i][2] *= scale;
    }
}
#endif

/**
 * @brief Draws the sample points to an OpenGL canvas.
 *
 * @param color The color to draw with.
 */
void MLSSurface::renderPoints(float *color) {

    glColor3fv(color);
    glPointSize(2.0);

    glBegin(GL_POINTS);
    for (int i = 0; i < m_VertexNr; i++) {
        
        if (m_InputNormals && (DIMENSION > 2))
            glNormal3f((GLfloat) m_Samples[i].n[0], (GLfloat) m_Samples[i].n[1], (GLfloat) m_Samples[i].n[2]);

        glVertex3f((GLfloat) m_Samples[i][0],(GLfloat)  m_Samples[i][1], (GLfloat) m_Samples[i][2]);
    }    
    glEnd();
    
}

/**
 * @brief Computes the MLS enery of a point given a vector field.
 *
 * @param query Point we want the energy for.
 * @param n Vector defined by the vector field.
 * @param kernel Defines support at the query point.
 *
 * Computes the enery of the query point.  This is defined as the
 * energy of a fitted plane to the weighted local support
 * (defined by kernel).
 */
double MLSSurface::computeMLSEnergy(const Point &query, Vector &n, Kernel *kernel)
{
    int i, found;
    double w,dist,totalweight=0;
    double energy=0;
    Point p; //temporary point so we don't have to call getPoint() all the time
    std::vector<DistancedIndex> neighbors;

    //The neighboring points as defined by the kernel
    kernel->getNeighbors(query, m_Tree, neighbors);
    
    found = neighbors.size();
    assert(found > 0); //We better have at least one neighbor out there somewhere
	if(found == 0) return 1.0e+100;

    //First we collect all weights and compute the weighted center of gravity
    for (i = 0; i < found; ++i) {
        p = m_Samples[neighbors[i].index];
        
        w = (*kernel)(query,p);
        totalweight += w;
        dist = dot(n, ((Point)p - query));
        energy += dist*dist*w;
    }   
    
    if (totalweight < FUDGE) {//If we are really far away we shouldn't trust the weighting
        //In this case we simply take our nearest neighbor (see the Amenta and Kil paper)
        std::nth_element(neighbors.begin(), neighbors.begin(), neighbors.end());
        p = m_Samples[neighbors.begin()->index];
        dist = dot(n, ((Point)p - query));
        return dist*dist;
    }
    
    return energy/totalweight;
}


/**
 * @brief Computes the Mahalanobis enery of a point given a vector field.
 *
 * @param query Point we want the energy for.
 * @param n Vector defined by the vector field.
 * @param kernel Defines support at the query point.
 *
 * TODO add more commenting
 */
double MLSSurface::computeMAHEnergy(const Point &query, Vector &n, Kernel *kernel)
{
    double w,dist,totalweight=0;
    double energy=0;
    int i, found;
    Point p; //temporary point so we don't have to call getPoint() all the time
    std::vector<DistancedIndex> neighbors;

    //The neighboring points as defined by the kernel
    kernel->getNeighbors(query, m_Tree, neighbors);    
    
    found = neighbors.size();
    assert(found > 0); //We better have at least one neighbor out there somewhere
	if(found == 0) return 1.0e+100;

    for (i = 0; i < found; ++i) {
        p = m_Samples[neighbors[i].index];
        
        w = (*kernel)(query,p);
        totalweight += w;
        dist = dot(n, ((Point)p - query));
        dist *= dist;
        dist += kernel->getMahFactor()*(((Point)p - query) - dot(n, ((Point)p - query))*n).lengthSquared();
        energy += dist*w;
    }   
    
    if (totalweight < FUDGE) {//If we are really far away we shouldn't trust the weighting
        //In this case we simply take our nearest neighbor (see the Amenta and Kil paper)
        std::nth_element(neighbors.begin(), neighbors.begin(), neighbors.end());
        p = m_Samples[neighbors.begin()->index];
        dist = dot(n, ((Point)p - query));
        dist *= dist;
        dist += kernel->getMahFactor()*(((Point)p - query) - dot(n, ((Point)p - query))*n).lengthSquared();
        return dist;
    }
    
    return energy/totalweight;
}
    
/**
 * @brief Computes the normal of a given point.
 *
 * @param query The point we wish to compute a normal for.
 * @param kernel Defines the support at query.
 * @param normal Return value: normal at query.
 * @param avg Center of gravity of the support.
 *
 * Computes the optimal direction using the Eigenvector corresponding to the larges Eigenvalue.
 */
bool MLSSurface::computeCovNormal(const Point &query,Kernel *kernel, Vector &normal, Point &avg)
{
    int i,j,found;
    double totalWeight = 0;
    Point p; //temporary point so we don't have to call getPoint() all the time
    std::vector<DistancedIndex> neighbors;

    //The neighboring points as defined by the kernel
    kernel->getNeighbors(query, m_Tree, neighbors);

    gsl_matrix *data_matrix;
    gsl_matrix *transposed_data_matrix;
    gsl_matrix *covariance_matrix;
    gsl_matrix *eigen_vectors;
    gsl_vector *eigen_values;
    gsl_eigen_symmv_workspace *workspace;
    
    found = neighbors.size();

    //If we found fewer points than necessary to set up the covariance matrix
    //this is pointless
    if (found < DIMENSION) {
        return false;
    }
    SCALARTYPE* weights = new SCALARTYPE[found]; // The set of weights

    avg = Point();  //sets it to 0,0,0: the origin.

    data_matrix = gsl_matrix_alloc(DIMENSION,found);
    transposed_data_matrix = gsl_matrix_alloc(found, DIMENSION);
    covariance_matrix = gsl_matrix_alloc(DIMENSION,DIMENSION);
    eigen_vectors = gsl_matrix_alloc(DIMENSION,DIMENSION);
    eigen_values = gsl_vector_alloc(DIMENSION);
    workspace = gsl_eigen_symmv_alloc(DIMENSION);
    
    //So now lets collect the (x - pi) values in the data matrix
    for (i=0;i<found;i++) {
        p = m_Samples[neighbors[i].index];
        
        for (j=0;j<DIMENSION;j++) 
            gsl_matrix_set(data_matrix,j,i,query[j] - p[j]);

        weights[i] = (float) (*kernel)(query,p);
        avg += weights[i]*(p);
        totalWeight += weights[i];
    }   
        
    for (j=0;j<DIMENSION;j++) 
        avg[j] /= totalWeight;
    
    // create transposed data matrix
    gsl_matrix_transpose_memcpy(transposed_data_matrix, data_matrix);
    
    // now we add the weighting factors to one matrix
    for (i=0;i<found;i++) {
        for (j=0;j<DIMENSION;j++) 
            gsl_matrix_set(data_matrix,j,i,gsl_matrix_get(data_matrix,j,i)*weights[i]);
    }
    
    delete[] weights;   // we're finished with weights

    // create matrix views
    gsl_matrix_const_view data_matrix_view = gsl_matrix_const_submatrix(data_matrix, 0,0, DIMENSION, found);
    gsl_matrix_const_view transposed_data_matrix_view = gsl_matrix_const_submatrix(transposed_data_matrix, 0,0, found, DIMENSION);
    // get covariance matrix
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &data_matrix_view.matrix, &transposed_data_matrix_view.matrix, 0.0, covariance_matrix);
    // calculate eigenvectors/eigenvalues of the covariance matrix
    
    gsl_matrix_set(covariance_matrix,0,0,gsl_matrix_get(covariance_matrix,0,0)/totalWeight);
    gsl_matrix_set(covariance_matrix,0,1,gsl_matrix_get(covariance_matrix,0,1)/totalWeight);
    gsl_matrix_set(covariance_matrix,1,0,gsl_matrix_get(covariance_matrix,1,0)/totalWeight);
    gsl_matrix_set(covariance_matrix,1,1,gsl_matrix_get(covariance_matrix,1,1)/totalWeight);
      
    gsl_eigen_symmv(covariance_matrix, eigen_values, eigen_vectors, workspace);
    gsl_eigen_symmv_sort(eigen_values, eigen_vectors, GSL_EIGEN_SORT_VAL_ASC);
    
    if (fabs(gsl_vector_get(eigen_values,0) - gsl_vector_get(eigen_values,1)) < 0.0001) {

        gsl_matrix_free(data_matrix);
        gsl_matrix_free(transposed_data_matrix);
        gsl_matrix_free(covariance_matrix);
        gsl_matrix_free(eigen_vectors);
        gsl_vector_free(eigen_values);
        gsl_eigen_symmv_free(workspace);
        
        return false;
    }
    
    // the normal vector is the one associated with the largest eigenvalue
    for (j=0;j<DIMENSION;j++) 
        normal[j] = gsl_matrix_get(eigen_vectors,j,0);
    normal.normalize();

    gsl_matrix_free(data_matrix);
    gsl_matrix_free(transposed_data_matrix);
    gsl_matrix_free(covariance_matrix);
    gsl_matrix_free(eigen_vectors);
    gsl_vector_free(eigen_values);
    gsl_eigen_symmv_free(workspace);

    return true;
}   
        
/**
 * @brief Computes the normal of a given point.
 *
 * @param query The point we wish to compute a normal for.
 * @param kernel Defines the support at query.
 * @param normal Return value: normal at query.
 * @param avg Center of gravity of the support.
 * @return True if successful.
 *
 * Computes the optimal direction from an MLS surface construction using the eigenvector corresponding
 * to the largest eigenvalue. The code is almost identical to "computeCOVNormal". However, here we
 * compute the covariance matrix as the sum of the difference not between query and the data point
 * but between the weighted average of query's neighbors and the data points. Essentially,
 * we compute the normal of the plane we would project query on.
 */
bool MLSSurface::computeOptMLSDir(const Point &query, Kernel *kernel, Vector &normal, Point &avg)
{
    int i,j,found;
    SCALARTYPE totalWeight = 0;
    Point p; //temporary point so we don't have to call getPoint() all the time
    std::vector<DistancedIndex> neighbors;

    //The neighboring points as defined by the kernel
    kernel->getNeighbors(query, m_Tree, neighbors);

    gsl_matrix *data_matrix;
    gsl_matrix *transposed_data_matrix;
    gsl_matrix *covariance_matrix;
    gsl_matrix *eigen_vectors;
    gsl_vector *eigen_values;
    gsl_eigen_symmv_workspace *workspace;
    
    found = neighbors.size();

    //If we found fewer points than necessary to set up the covariance matrix
    //this is pointless
    if (found < DIMENSION) {
        return false;
    }

    SCALARTYPE* weights = new SCALARTYPE[found]; // The set of weights
    avg = Point();  //sets it to 0,0,0: the origin.

    data_matrix = gsl_matrix_alloc(DIMENSION,found);
    transposed_data_matrix = gsl_matrix_alloc(found, DIMENSION);
    covariance_matrix = gsl_matrix_alloc(DIMENSION,DIMENSION);
    eigen_vectors = gsl_matrix_alloc(DIMENSION,DIMENSION);
    eigen_values = gsl_vector_alloc(DIMENSION);
    workspace = gsl_eigen_symmv_alloc(DIMENSION);
    
    //First we collect all weights and compute the weighted center of gravity
    for (i=0;i<found;i++) {
        p = m_Samples[neighbors[i].index];
        
        weights[i] = (float) (*kernel)(query,p);
        avg += weights[i]*(p);
        totalWeight += weights[i];
    }   
    
    for (j=0;j<DIMENSION;j++) 
        avg[j] /= totalWeight;
    
    //So now lets collect the (x - pi) values in the data matrix
    for (i=0;i<found;i++) {
        p = m_Samples[neighbors[i].index];
        
        for (j=0;j<DIMENSION;j++) 
            gsl_matrix_set(data_matrix,j,i,(avg[j] - p[j]));
    }   
        
    // create transposed data matrix
    gsl_matrix_transpose_memcpy(transposed_data_matrix, data_matrix);
    
    
    // now we add the weighting factors to one matrix
    for (i=0;i<found;i++) {
        for (j=0;j<DIMENSION;j++) 
            gsl_matrix_set(data_matrix,j,i,gsl_matrix_get(data_matrix,j,i)*weights[i]);
    }

    delete[] weights;   // we're finished with weights
    
    // create matrix views
    gsl_matrix_const_view data_matrix_view = gsl_matrix_const_submatrix(data_matrix, 0,0, DIMENSION, found);
    gsl_matrix_const_view transposed_data_matrix_view = gsl_matrix_const_submatrix(transposed_data_matrix, 0,0, found, DIMENSION);
    // get covariance matrix
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &data_matrix_view.matrix, &transposed_data_matrix_view.matrix, 0.0, covariance_matrix);
    // calculate eigenvectors/eigenvalues of the covariance matrix
    
    gsl_matrix_set(covariance_matrix,0,0,gsl_matrix_get(covariance_matrix,0,0)/totalWeight);
    gsl_matrix_set(covariance_matrix,0,1,gsl_matrix_get(covariance_matrix,0,1)/totalWeight);
    gsl_matrix_set(covariance_matrix,1,0,gsl_matrix_get(covariance_matrix,1,0)/totalWeight);
    gsl_matrix_set(covariance_matrix,1,1,gsl_matrix_get(covariance_matrix,1,1)/totalWeight);
    
    gsl_eigen_symmv(covariance_matrix, eigen_values, eigen_vectors, workspace);
    gsl_eigen_symmv_sort(eigen_values, eigen_vectors, GSL_EIGEN_SORT_VAL_ASC);
    
    if (fabs(gsl_vector_get(eigen_values,0) - gsl_vector_get(eigen_values,1)) < 0.0001) {
        gsl_matrix_free(data_matrix);
        gsl_matrix_free(transposed_data_matrix);
        gsl_matrix_free(covariance_matrix);
        gsl_matrix_free(eigen_vectors);
        gsl_vector_free(eigen_values);
        gsl_eigen_symmv_free(workspace);
        
        return false;
    }

    // the normal vector is the one associated with the largest eigenvalue
    for (j=0;j<DIMENSION;j++) 
        normal[j] = gsl_matrix_get(eigen_vectors,j,0);
    normal.normalize();

    gsl_matrix_free(data_matrix);
    gsl_matrix_free(transposed_data_matrix);
    gsl_matrix_free(covariance_matrix);
    gsl_matrix_free(eigen_vectors);
    gsl_vector_free(eigen_values);
    gsl_eigen_symmv_free(workspace);

    return true;
}   
        
/**
 * @brief Computes the normal of a given point.
 *
 * @param query The point we wish to compute a normal for.
 * @param kernel Defines the support at query.
 * @param normal Return value: normal at query.
 * @param avg Center of gravity of the support.
 * @return True if successful.
 *
 * Computes the normal at the query simply by averaging the weighted normals
 * of the surflet support defined by kernel.
 */
bool MLSSurface::computeAvgNormal(const Point &query, Kernel *kernel, 
    Vector &normal, Point &avg)
{
    int i,found;
    Point p;
    SCALARTYPE totalWeight = 0,w;
    normal = Vector();  // make a zero vector
    std::vector<DistancedIndex> neighbors;

    //The neighboring points as defined by the kernel
    kernel->getNeighbors(query, m_Tree, neighbors);

    found = neighbors.size();

    //If we found no neighbor this is pointless
    if (found <= 0)
        return false;
    
    //Now collect all weighted normals
    for (i=0;i<found;i++) {
        p = m_Samples[neighbors[i].index];
        w = (float) (*kernel)(query,p);
        avg += w*(p);
        normal += w*m_Samples[neighbors[i].index].n;
        totalWeight += w;
    }

    normal *= 1/totalWeight;
    normal.normalize();
    
    return true;
}

bool MLSSurface::computeDistanceAvg(const Point &query, Kernel *kernel, double &distAvg)
{
    int found;
    std::vector<DistancedIndex> neighbors;

    //The neighboring points as defined by the kernel
    kernel->getNeighbors(query, m_Tree, neighbors);
        
    found = neighbors.size();

    //If we found no neighbor this is pointless
    if (found <= 0)
        return false;
    
    distAvg = 0;
    SCALARTYPE totWeight = 0, w;
    Sample s;
    //compute the average distance to the sample points
    for (int i = 0; i < found; ++i) {
        s = m_Samples[neighbors[i].index];
        w = (float) (*kernel)(query,s);

        distAvg += w * dot(s.n, (query - s));
        totWeight += w;
    }   

    distAvg /= totWeight; //This calculation normalizes the distance
    return true;
}

/**
 * @brief Projects a point onto the surface.
 *
 * @param query The point to project.
 * @param kernel Defines the support for the point.
 * @param optDir A function defining a vector field.
 * @param energy A function defining an energy for a point given a normal.
 * @param proj Return value. The projected location on the surface.
 * @return True if the projection succeeded.
 *
 * This function projects a point onto the surface as described
 * in "Defining Point-Set Surfaces" by Amenta and Kil.
 */
bool MLSSurface::projectPoint(const Point &query, Kernel *kernel, 
    bool (MLSSurface::*optDir)(const Point&, Kernel*, Vector&, Point&),
    double (MLSSurface::*energy)( const Point&, Vector&, Kernel*),
    Point &proj)
{                                                               
    Vector normal;
    Point avg;
    Point current,last;
    bool flag;
    int iter = 0;
    SCALARTYPE t;

    current = query;
        
    do {
        iter++;
        ///First we compute the optimal direction
        flag = ((*this).*optDir)(current,kernel,normal,avg);
    
        /// if we cannot compute a direction we bail out
        if (flag == false) {
            proj = current;
            return false;
        }
        t = (float) this->findEnergyMin(current,kernel,normal,MAX(sqrt(kernel->getKernelSizeSquared()), distance(avg,current)),energy);
        last = current;
        current = current + t*normal;

    } while ((distanceSquared(current,last) > SURFACE_TOL) && (iter < PROJ_STEPS));

    /*if (iter >= PROJ_STEPS) 
        {printf("failed to converge\n"); fflush(stdout);}*/

    //NB: returning false means we did not converge.  However, we still return some estimate
    //of the projected point.
    proj = current;
    return iter < PROJ_STEPS;
}

/**
 * @brief Finds the minimum of the an energy function.
 *
 * @param query Base point for the plane defining the energy.
 * @param kernel Defines the support around the query.
 * @param normal Normal of the plane defining the energy.
 * @param bound Bounds for the test.
 * @param energy Energy function that we will minimize.
 * @return Offset from query along the normal that is the first local minima
 *         of the function energy.
 *
 * Computes the first local minima of the function energy.  The domain of our search
 * is the line defined by normal.  We begin at query and search for the first minima
 * we encounter.  When we find it we return this as an offset along the normal.
 */
double MLSSurface::findEnergyMin(const Point &query, Kernel *kernel, Vector &normal,
    double bound, double (MLSSurface::*energy)(const Point&, Vector&, Kernel*))
{
    
    ///Setup the local variables to fit the gsl interface
    static ParamStruct minStruct;
    minStruct.surface = this;
    minStruct.f = energy;
    minStruct.q = query;
    minStruct.n = normal;
    minStruct.k = kernel;

    int status, iter=0, max_iter = 20;
    double left=-bound,right=bound,min = 0; ///boundaries and initial guess of the minimization
    const gsl_min_fminimizer_type *T;
    gsl_min_fminimizer *s;
    gsl_function F; ///gsl internal format
    double l,m,r;
    
    l = gslMinFuncInterface(-bound,&minStruct);
    m = gslMinFuncInterface(0,&minStruct);
    r = gslMinFuncInterface(bound,&minStruct);
    
    //In case we are far away from the surface "bound" might not be wide enough
    //so as a security we than go into the direction of smallest energy
    if ((l <= m) && (l <= r))
        return -bound;
    else if ((r <= m) && (r < l))
        return bound;
    
    //return 0;
    F.function = &gslMinFuncInterface; 
    F.params = &minStruct;

    T = gsl_min_fminimizer_brent;
    s = gsl_min_fminimizer_alloc (T);
    gsl_min_fminimizer_set (s, &F, min, left, right);

    do {
        iter++;
        status = gsl_min_fminimizer_iterate (s);

        min = gsl_min_fminimizer_x_minimum (s);
        left = gsl_min_fminimizer_x_lower (s);
        right = gsl_min_fminimizer_x_upper (s);
        
        status = gsl_min_test_interval (left, right, 0.0001, 0.0);
    } while (status == GSL_CONTINUE && iter < max_iter);
    
    gsl_min_fminimizer_free(s);
    return 0.5*(left+right);
}

/**
 * @brief Returns the MLS's kernel.
 */
Kernel *MLSSurface::getKernel() {

    return m_Kernel;
}

/**
 * @brief Function we want GSL to minimize.
 */
double MLS::gslMinFuncInterface(double x, void *param) { 

    ParamStruct *min = (ParamStruct*)param;
    Point current;
    current = min->q + ((SCALARTYPE)x)*min->n;
    
    return ((*min->surface).*(min->f))(current,min->n,min->k);
}


/**
 * @brief Kernel Selector combobox implementation.
 */

MLSSurface::KernelSelector::KernelSelector(Surface * surface, Kernel** f)
:SurfParamComboBox::Callback()
,_basicKernel(f)
{
	//We could replace this by an own factory... but to me it seems overkill...
	_choices.push_back(std::string("Gauss"));
	_choices.push_back(std::string("Gaussian center and linear tails to 0"));
	(*_basicKernel)= new GaussKernel();
}

MLSSurface::KernelSelector::~KernelSelector()
{
	delete (*_basicKernel);
	(*_basicKernel)=0;
}

void MLSSurface::KernelSelector::itemselected(unsigned int selection)
{
	if (selection==0)
	{
		delete (*_basicKernel);
		(*_basicKernel)=new GaussKernel();
	}
	else
	{
		delete (*_basicKernel);
		(*_basicKernel)=new GaussianClip();		
	}
}

#endif