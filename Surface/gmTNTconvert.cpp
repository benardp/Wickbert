/**
 * This file contains several global "convert" methods to easily convert 
 * gmVector3/4 <-> TNT::Vector and gmMatrix3/4 <-> TNT::Matrix.  
 * @file gmTNTconvert.cpp
 * @author Terry Fleury <tfleury@uiuc.edu>
 * @date October 28, 2001
 */

#include "gmTNTconvert.h"

/**
 * Convert a gmVector3 to a TNT::Vector<double> of size 3.
 * @param v The input gmVector3 to convert.
 * @return The equivalent TNT::Vector.
 */
TNT::Vector<double> convert(const gmVector3 & v)
{
  TNT::Vector<double> tv(3);
  for (int i = 0; i < 3; i++)
    tv[i] = v[i];
  return tv;
}

/**
 * Convert a gmMatrix3 to a TNT::Matrix<double> of size 3x3.
 * @param m The input gmMatrix3 to convert.
 * @return The equivalent TNT::Matrix.
 */
TNT::Matrix<double> convert(gmMatrix3 m)
{
  TNT::Matrix<double> tm(3,3);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      tm[i][j] = m[i][j];
  return tm;
}

/**
 * Convert a gmVector4 to a TNT::Vector<double> of size 4.
 * @param v The input gmVector4 to convert.
 * @return The equivalent TNT::Vector.
 */
TNT::Vector<double> convert(gmVector4 v)
{
  TNT::Vector<double> tv(4);
  for (int i = 0; i < 4; i++)
    tv[i] = v[i];
  return tv;
}

/**
 * Convert a gmMatrix4 to a TNT::Matrix<double> of size 4x4.
 * @param m The input gmMatrix3 to convert.
 * @return The equivalent TNT::Matrix.
 */
TNT::Matrix<double> convert(gmMatrix4 m)
{
  TNT::Matrix<double> tm(4,4);
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      tm[i][j] = m[i][j];
  return tm;
}

/** 
 * Convert a TNT::Vector<double> of size 3 to a gmVector3.
 * @param v The input TNT::Vector to convert.
 * @return The equivalent gmVector3.
 * @note If the TNT::Vector is not of size 3, then only the available
 * elements of the TNT::Vector are copied, thus avoiding invalid array
 * accesses.
 */
gmVector3 convert(TNT::Vector<double> v)
{
  gmVector3 gv;
  for (int i = 0; i < v.size(); i++)
    gv[i] = v[i];
  return gv;
}

/** 
 * Convert a TNT::Matrix<double> of size 3x3 to a gmMatrix3.
 * @param m The input TNT::Matrix to convert.
 * @return The equivalen gmMatrix3.
 * @note If the TNT::Matrix is not of size 3x3, then only the available
 * elements of the TNT::Matrix are copied, thus avoiding invalid array
 * accesses.
 */
gmMatrix3 convert(TNT::Matrix<double> m)
{
  gmMatrix3 gm;
  for (int i = 0; i < m.num_rows(); i++)
    for (int j = 0; j < m.num_cols(); j++)
      gm[i][j] = m[i][j];
  return gm;
}

/** 
 * Convert a TNT::Vector<double> of size 4 to a gmVector4.
 * @param v The input TNT::Vector to convert.
 * @return The equivalent gmVector4.
 * @note If the TNT::Vector is not of size 4, then only the available
 * elements of the TNT::Vector are copied, thus avoiding invalid array
 * accesses.
 * @note This method is named convert4 to disambiguate it from
 * convert(TNT::Vector<double> &v), which returns a gmVector3 instead of a
 * gmVector4.  Unfortunately, there's no way to name them both "convert"
 * since they have the same arguments.
 */
gmVector4 convert4(TNT::Vector<double> v)
{
  gmVector4 gv;
  for (int i = 0; i < v.size(); i++)
    gv[i] = v[i];
  return gv;
}

/** 
 * Convert a TNT::Matrix<double> of size 4x4 to a gmMatrix4.
 * @param m The input TNT::Matrix to convert.
 * @return The equivalent gmMatrix4.
 * @note If the TNT::Matrix is not of size 4x4, then only the available
 * elements of the TNT::Matrix are copied, thus avoiding invalid array
 * accesses.
 * @note This method is named convert4 to disambiguate it from
 * convert(TNT::Matrix<double> &m), which returns a gmMatrix3 instead of a
 * gmMatrix4.  Unfortunately, there's no way to name them both "convert"
 * since they have the same arguments.
 */
gmMatrix4 convert4(TNT::Matrix<double> m)
{
  gmMatrix4 gm;
  for (int i = 0; i < m.num_rows(); i++)
    for (int j = 0; j < m.num_cols(); j++)
      gm[i][j] = m[i][j];
  return gm;
}

