/**
 * This file contains the method definitions for a CriticalMesh object.
 * @file   CriticalMesh.cpp
 * @date   April 22, 2002
 * @author Terry Fleury (tfleury@uiuc.edu)
 */

#include "CriticalMesh.h"

void CriticalMesh::initialMesh(Box3d bounds)
{
  if (getSurface() == NULL)
    return;  // Do nothing if we don't have an Implicit to work with

  // STEP 1 - Find and classify all critical points in the bounding box.
  setBoxBounds(bounds);  
  fUseNewton = true;
  search();
  classify();

  // STEP 2 - 
}

