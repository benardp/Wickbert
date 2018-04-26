/*!
 * Declaration of an abstract displayable/savable object
 * \file GlobalListObject.h
 * \date September 19, 2001
 * \author Wojciech Jarosz
 */

#ifndef GLOBALLISTOBJECT_H
#define GLOBALLISTOBJECT_H

#include "Implicit.h"

/*!
 * A Global List Object stores a pointer to an Implicit Object,
 * a display flag, and a save flag.
 */
class GlobalListObject
{
  public:
    Implicit* surface; //!< Pointer to the Implicit surface
    bool display;      //!< True/False flag should object be displayed somehow 
    bool save;         //!< True/False flag should object be saved

    /// Constructor
    GlobalListObject(Implicit* mysurface = NULL, 
                     const bool& mydisplay = true, 
                     const bool& mysave = true)
    {
      surface  = mysurface;
      display  = mydisplay;
      save  = mysave;
    }

    /// Destructor
    ~GlobalListObject()
    {
      delete surface;
    }
};

#endif

