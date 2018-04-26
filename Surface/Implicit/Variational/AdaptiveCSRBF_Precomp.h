/** @file AdaptiveCSRBF_Precomp.h
 * Same as AdaptiveCSRBF, except loaded file is expected to be a precomputed set of approximation
 * centers. Cannot be modified at runtime (no parameters)
 * @author: Scott Kircher
 * @date: December 14, 2004
 */

#ifndef ADAPTIVECSRBF_PRECOMP_H
#define ADAPTIVECSRBF_PRECOMP_H

#include "Surface/Implicit/Variational/AdaptiveCSRBF.h"

/** AdaptiveCSRBF_Precomp class
 */
class AdaptiveCSRBF_Precomp : public AdaptiveCSRBF
{
public:
	MAKE_NAME();

	// Does nothing
	virtual void updateRBF(void) {;}

	unsigned int qlen() {return 0;}
	void getq(double *q) {;}
	void _setq(double *q) {;}
	void getqname(char **qn) {;}

	/// read the precomputed approximation centers in from a file
	virtual bool readImplicit(std::ifstream &file, bool verbose);
};

#endif