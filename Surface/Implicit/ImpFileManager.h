/***
 * Declaration of the implicit file manager class.
 * @file ImpFileManager.h
 * @date July 11, 2001
 * @author Ed Bachta
 */

#ifndef __ASL_IMPFILEMANAGER_H__
#define __ASL_IMPFILEMANAGER_H__

#include "Implicit.h"

#include <fstream>
#include <iostream>
#include <string>

/**
 * The file manager class allows for reading and writing implicit 
 * surface descriptions to *.imp files. It makes use of factories to 
 * create objects from file descriptions.
 *
 * The file manager has a version number and checks this against the
 * version number that appears at the top of .imp files. If the file
 * format ever changes this can be used to read the file properly based
 * on version number.
 *
 * Usage:
 *
 * ImpFileManager ifm; // your very own file manager
 * Implicit* surface;  // your handy surface pointer
 *
 * ifm.readFile("surf.imp", surface);  // read from file
 * ifm.writeFile("surf.imp", surface); // write to file
 *
 * You can also leave off the .imp extension and the file manager will
 * add it in automatically.
 */
class ImpFileManager {

 private:

  const std::string suffix;  ///< Implicit file suffix.
  const std::string version; ///< Version of this file manager.
	
  int verbose;

 public:

  /**
   * Constructor.
   */
  ImpFileManager();

  /**
   * Set verbose mode.
   * @li 0 - no output
   * @li 1 - standard output
   * @li 2 - light debug output
   */
  void setVerbose(int v);

  /// Read an implicit object from file.
  bool readImplicit(std::ifstream &file,
                    Implicit * &object,
                    std::string &name);

  /// Write an implicit object to file.
  bool writeImplicit(std::ofstream &file,
                     Implicit * &object,
                     std::string tabs);

  /// Read an implicit surface from file.
  bool readFile(const std::string& name, Implicit* &surface);

  /// Read a group of implicit surfaces from file.
  bool readFile(const std::string& name, std::vector<Surface*>& surfaces);

  /// Write a group of implicit surfaces to file.
  bool writeFile(const std::string& name, std::vector<Surface*>& surfaces);

  /// Write an implicit surface to file.
  bool writeFile(const std::string& name, Implicit* &surface);

};

#endif

