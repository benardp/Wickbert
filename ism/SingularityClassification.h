/**
 * Declaration of a class used to classify the critical points of a feature contour.
 * @file SingularityClassification.h
 * @date 10/14/2005
 * @author Matei N. Stroila
 * @remarks See Simon Plantinga, Gert Vegter "Contour generators of evolving implicit surfaces" 
 * http://doi.acm.org/10.1145/781606.781614 .
 */

#include "Newton.h"

class Implicit;

class SingularityClassification{
	
public:
	
	/**
	 * Constructor.
	 * @param rotF rotated implicit (such that the view vector is (0,0,1))
	 */
	SingularityClassification(Implicit *rotF);
	~SingularityClassification(){};
	
	/**
	 * Function that returns the singularity type of a critical point positioned at myP.
	 * @param myP the postion vector
	 */
	SingularityType getSingularityType(gmVector3& myP);
	
private:
	
	///The sign of sigma determines the singularity type
	double Sigma; 
	///Rotated implicit (such that the view vector is (0,0,1)). 
	Implicit* F; 
	///Position vector of a singular point of the contour.
	gmVector3 p; 

	void computeSigma();
	
};
