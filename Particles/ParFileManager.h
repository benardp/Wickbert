/***
 * Declaration of the particle file manager class.
 * @file ParFileManager.h
 * @date July 6, 2002
 * @author John C. Hart
 * @note Based on ImpFileManager.h
 */

#ifndef __ASL_PARFILEMANAGER_H__
#define __ASL_PARFILEMANAGER_H__

#include "Particles.h"

#include <fstream>
#include <iostream>
#include <string>

/**
 * The file manager class allows for reading and writing particle systems
 * to *.par files. It makes use of factories to 
 * create objects from file descriptions.
 *
 * The file manager has a version number and checks this against the
 * version number that appears at the top of .par files. If the file
 * format ever changes this can be used to read the file properly based
 * on version number.
 *
 * Usage:
 *
 * ParFileManager pfm; // your very own file manager
 * std::vector<Particles*> plist;  // your handy list of particle systems
 *
 * ifm.readFile("surf.imp", plist);  // read from file
 * ifm.writeFile("surf.imp", plist); // write to file
 *
 * You can also leave off the .par extension and the file manager will
 * add it in automatically.
 */
class ParFileManager {

 private:

  const std::string suffix;  ///< Implicit file suffix.
  const std::string version; ///< Version of this file manager.
	
  int verbose;

protected:
  virtual bool pfmerror(std::string s) {
    std::cerr << s << std::endl;
    return false;
  }

  virtual void pfmstatus(int v, std::string s) {
    if (verbose >= v)
      std::cerr << s << std::endl;
  }



 public:

  /**
   * Constructor.
   */
  ParFileManager();

  ///virtual destructor to avoid the annoying warnings
  virtual ~ParFileManager() {}
	
  /**
   * Set verbose mode.
   * @li 0 - no output
   * @li 1 - standard output
   * @li 2 - light debug output
   */
  void setVerbose(int v);

  /// Read a Particles object from file.
  bool readParticles(std::ifstream &file, Particles *&p, const std::string& name);

  /// Write a Particles object to file.
  bool writeParticles(std::ofstream &file, Particles *p);

  /// Read an implicit surface from file.
  bool readFile(const std::string& name, Particles* &p);

  /// Write an implicit surface to file.
  bool writeFile(const std::string& name, Particles* &p);

  /// Read a group of implicit surfaces from file.
  bool readFile(const std::string& name, std::vector<Particles*>& p);

  /// Write a group of implicit surfaces to file.
  bool writeFile(const std::string& name, std::vector<Particles*>& p);

  /// write par particle stuff to an obj file
  bool writePerParticle(const std::string& name, Particles *p);

  /// called by previous, this is a helper
  bool writePerParticle(std::ofstream &file, Particles *p);
};

#endif
