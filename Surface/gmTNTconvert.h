/**
 * This file contains function prototypes for several global "convert"
 * methods to easily convert gmVector3/4 <-> TNT::Vector and 
 * gmMatrix3/4 <-> TNT::Matrix.  
 * @file gmTNTconvert.h
 * @author Terry Fleury <tfleury@uiuc.edu>
 * @date October 28, 2001
 */

#ifndef GM_TNT_CONVERT_H
#define GM_TNT_CONVERT_H

#include "tnt/cmat.h"
#include "tnt/vec.h"
#include "libgm/gm.h"

// Convert a gmVector3 to a TNT::Vector<double> of size 3.
TNT::Vector<double> convert(const gmVector3 & );
// Convert a gmMatrix3 to a TNT::Matrix<double> of size 3x3.
TNT::Matrix<double> convert(gmMatrix3);
// Convert a gmVector4 to a TNT::Vector<double> of size 4.
TNT::Vector<double> convert(gmVector4);
// Convert a gmMatrix4 to a TNT::Matrix<double> of size 4x4.
TNT::Matrix<double> convert(gmMatrix4);
// Convert a TNT::Vector<double> of size 3 to a gmVector3.
gmVector3 convert(TNT::Vector<double>);
// Convert a TNT::Matrix<double> of size 3x3 to a gmMatrix3.
gmMatrix3 convert(TNT::Matrix<double>);
// Convert a TNT::Vector<double> of size 4 to a gmVector4.
gmVector4 convert4(TNT::Vector<double>);
// Convert a TNT::Matrix<double> of size 4x4 to a gmMatrix4.
gmMatrix4 convert4(TNT::Matrix<double>);

#endif

