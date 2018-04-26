#include "Surface.h"
#include <iostream>

/** Computes the outer product between to 3-vectors. The result might
 * be the transpose of what it really should be. I haven't checked the
 * ordering.
 * @param x 3-vector factor
 * @param y 3-vector factor
 * @return 3x3 matrix containing pairwise products of x,y elements
 * @note Hate to put this here of all places, but want to keep libgm unchanged
 * and this is the only file I know will be included by everyone.
 */
gmMatrix3 outer( gmVector3 x, gmVector3 y)
{
  gmMatrix3 xy(x[0]*y[0],x[1]*y[0],x[2]*y[0],
               x[0]*y[1],x[1]*y[1],x[2]*y[1],
               x[0]*y[2],x[1]*y[2],x[2]*y[2]);
  return xy;
}

/// Interval version of outer product of two Boxes
IMatrix3d outer(Box3d x, Box3d y)
{
  IMatrix3d xy(x[0]*y[0],x[1]*y[0],x[2]*y[0],
               x[0]*y[1],x[1]*y[1],x[2]*y[1],
               x[0]*y[2],x[1]*y[2],x[2]*y[2]);
  return xy;
}



/// default constructor
Surface::Surface()
:objectName("")
{}

// Definition of static member string defaultObjectNamePrefix
std::string Surface::defaultObjectNamePrefix = std::string("~!");
unsigned int Surface::objectNumber = 0;

/**
 * Return the name of this Surface.  Since Surface does not define
 * registry_name in Surface.cpp, it returns NULL.  This should be the
 * case for all abstract classes.
 * @return registry_name The name of the Surface class.
 * @todo Make sure that NULL is the default value of registry_name in
 * Unix/Linux.
 * @todo This member is simply "std::string name" for Particles. Should
 * remove the function and just return the string.
 */
const std::string Surface::name() 
{
  return registry_name; 
}

/**
 * Returns the name of the instanced object WITHOUT any 'default' prefix. 
 * If there is no 'default' prefix, then this method returns the same string
 * as getFullObjectName(). This method returns a std::string.
 * @return The name of the Surface object as a string sans default prefix.
 * @see getFullObjectName();
 * @see getObjectNameAsCString();
 */
std::string Surface::getObjectName()
{
  checkObjectName();
  if (objectNameIsDefault())
    return objectName.substr(defaultObjectNamePrefix.length());
  else
    return objectName;
}

/**
 * Returns the name of the object (sans prefix) as a C-string.  If there is
 * no 'default' prefix, then this method returns the same string as
 * getFullObjectNameAsCString().  This method returns a C-string (char*).
 * @return The name of the Surface object as a C-string sans default prefix.
 * @see getFullObjectNameAsCString();
 * @see getObjectName();
 */
const char* Surface::getObjectNameAsCString()
{
  return getObjectName().c_str();
}

/** 
 * Returns the name of the instanced object WITH any 'default' prefix.  If
 * there is no 'default' prefix, then this method returns the same string as
 * getObjectName().  This method returns a std::string.
 * @return The name of the Surface object as a String with default prefix.
 * @see getObjectName();
 * @see getFullObjectNameAsCString();
 */
std::string Surface::getFullObjectName()
{
  checkObjectName();
  return objectName;
}

/**
 * Returns the name of the object (with prefix) as a C-string.  If there is
 * no 'default' prefix, then this method returns the same string as
 * getObjectNameAsCString().  This method returns a C-string (char*).
 * @return The name of the Surface object as a C-string with default prefix.
 * @see getObjectNameAsCString();
 * @see getFullObjectName();
 */
const char* Surface::getFullObjectNameAsCString()
{
  return getFullObjectName().c_str();
}

/**
 * Sets the name of the object to the passed-in string.  Note that you
 * CANNOT set the name of the object to a 'default' state by passing in a
 * string with the 'default' prefix.  This method will remove any 'default'
 * prefix.  If you want a default name for the object, use
 * resetObjectName() OR pass in an empty string (see Note). 
 * @param newName The new name for the object.
 * @see resetObjectName();
 * @note If you pass in an empty string, this method will set the string to
 * be a 'default' name without any value. This will erase any name currently
 * assigned to the object AND make the string equal to the 'default' prefix
 * so that you don't save the empty name string to file.
 * @param newName The new name for the object.
 * @see setObjectName(char*);
 */
void Surface::setObjectName(std::string newName)
{
  if (newName.empty())
    objectName = defaultObjectNamePrefix;
  else 
    {
      objectName = newName;
      // Remove any default prefixes from the objectName
      while (objectNameIsDefault())
        {
          objectName.erase(0,defaultObjectNamePrefix.length());
        }
    }
}

/**
 * Sets the name of the object to the passed-in C-string.  This method is
 * the same as setObjectName(std::string) but takes a C-string as the
 * parameter.
 * @param newName The new name for the object.
 * @see setObjectName(std::string)
 */
void Surface::setObjectName(char* newName)
{
  setObjectName(std::string(newName));
}
   
/// Sets the name of the object to a unique 'default' name.
void Surface::resetObjectName() 
{ 
  objectName = defaultObjectNamePrefix; 
}

/**
 * Returns true if the object name has a 'default' prefix.  This method is
 * useful in determining if the object's name is in a 'default' format so
 * that you don't save the object's name to file.
 * @return True if the object's name is in a default format (ie. has a
 *         'default' prefix).
 */
bool Surface::objectNameIsDefault()
{
  return (objectName.substr(0,defaultObjectNamePrefix.length()) == 
    defaultObjectNamePrefix);
}

/**
 * Sets the internal unique object counter number.  This should only be used
 * if you know what you are doing because it affects the names of all newly
 * created Surface objects.  You may want to use this after you have read
 * in an Surface object file and notice that somehow one of the names is in
 * the default object name format and need to update the internal object
 * counter to be higher than that object's number.  Use with care!
 * @param newNumber The new number for the internal object counter.
 */
void Surface::setObjectNumber(int newNumber)
{
  if (newNumber >= 0)
    objectNumber = newNumber;
}

/**
 * Checks to see if the objectName has not been set yet.  If the objectName
 * value is empty (ie. ""), then we call resetObjectName to give the
 * object a 'default' object name.  This method is called every time you 
 * try to get the objectName to make sure that the objectName has been
 * properly initialized.  The objectName string can never be empty since
 * if you try to pass the empty string to setObjectName(), that method
 * resets the object name to be the default prefix.  If the objectName
 * string is ever empty, then something has gone wrong!
 */
void Surface::checkObjectName()
{
  if (objectName.empty())
    resetObjectName();
}

std::ostream &operator<<(std::ostream &out, Surface *s)
{
	std::string nm = s->name(); // does this name exist in the class factory?
	std::string gfon = s->getFullObjectName();

	// output class name
	out << s->name() << " ";

	// output name of this instance
	out << '"' << s->getFullObjectName() << "\" {\n";

	out << &(s->params);

	/* output type specific info (e.g. q's in Implicit)
	 * using virtual output() routine */
	s->output(out);

	out << "}\n\n" << std::endl;

	return out;
}

std::ostream &operator<<(std::ostream &out, Surfaces *s)
{
	for (unsigned int i = 0; i < s->size(); i++)
		out << s->at(i);

	return out;
}

