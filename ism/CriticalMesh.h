/**
 * This file contains the class declaration for a CriticalMesh object.
 * @file   CriticalMesh.h
 * @date   April 22, 2002
 * @author Terry Fleury (tfleury@uiuc.edu)
 */

#ifndef CRITICALMESH_H
#define CRITICALMESH_H

#include "SearchCritical.h"

/**
 * 
 */
class CriticalMesh : public SearchCritical 
{
  public:
    void initialMesh(Box3d);
};

#endif

