/**
 * @file ParticleEigen.h
 * @author Terry Fleury <tfleury@uiuc.edu>
 * This attribute stores the eigenvectors, eigenvalues, and number of
 * positive eigenvalues of a point. This is useful in classifying
 * the type of critical point.
 * 
 * The type of critical point is based on the number of negative eigenvalues
 * at each point.  A MAX critical point has three decreasing eigen vectors
 * ( = 3 negative eigen values ) and a MIN critical point has three
 * increasing eigen vectors ( = 3 positive eigen values ).  Note that these
 * are labelled 'backwards' from John C. Hart's presentation on critical
 * points.
 * Positive eigenvalues     Critical Point Type
 * --------------------     -------------------
 *          0                    Maximum
 *          1                 Saddle Type 2
 *          2                 Saddle Type 1
 *          3                    Minimum
 */

#ifndef PARTICLEEIGEN_H
#define PARTICLEEIGEN_H

#include "Particles.h"
#include "ParticleAttribute.h"

class ParticleEigen : public ParticleAttribute 
{
protected:
	/// Eigen vector
  std::vector<gmMatrix3> eigenVectors;
  /// Eigen value
  std::vector<gmVector3> eigenValues;
  /// Number of positive eigen values => critical point type
  std::vector<int>positiveEigenValues;

public:
	MAKE_PARTICLESTUFF_NAME();

	/// Add paritcle orientation to a system of particles.
	ParticleEigen(Particles *ps=NULL, const std::string& name=std::string("ParticleEigen"));

	/// Callback for particle addition.
	void particleAdded();

	/// Callback for particle removal.
	void particleRemoved(unsigned int i);

	/// Getters/setters for the protected class variables
	virtual gmMatrix3 getEigenVectors(unsigned int);
	virtual void setEigenVectors(unsigned int, gmMatrix3);
  virtual gmVector3 getEigenValues(unsigned int);
  virtual void setEigenValues(unsigned int, gmVector3);
  virtual int getPositiveEigenValues(unsigned int);
  virtual void setPositiveEigenValues(unsigned int, int);

  virtual void clear();
};

#endif

