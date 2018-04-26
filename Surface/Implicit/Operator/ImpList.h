/**
 * Declaration of the implicit list.
 * @file ImpList.h
 * @date July 18, 2001
 * @author Ed Bachta
 */

#ifndef __ASL_IMPLIST_H__
#define __ASL_IMPLIST_H__

#include "Surface/Implicit/Implicit.h"

#include <list>

/**
 * The ImpList class is a list of implicit surfaces. f(x) for
 * an implicit list is defined here as the summation over the list
 * of each implicit's f(x).
 */
class ImpList : public Implicit, public std::list<Implicit*>
{
  private:
    std::list<Implicit*>::iterator i_it;   ///< Convenience iterator variable

  public:
    ImpList();                             ///< Default constructor.
    ImpList(std::list<Implicit*> ilist);   ///< Explicit constructor.

#ifndef INTERVAL_EVAL_ONLY
    virtual double proc(const gmVector3 & );        ///< Evaluate the sum of f(x).
    virtual gmVector3 grad(const gmVector3 & );     ///< Evaulate the sum of grad(x).
    virtual gmMatrix3 hess(const gmVector3 & );     ///< Evaluate the sum of hess(x).
#endif

    virtual Intervald proc(const Box<double>& );
    virtual Box3d grad(const Box<double>& );
    virtual IMatrix3d hess(const Box<double>& );

    virtual unsigned int qlen();                    ///< Return number of parameters.
    virtual void getq(double*);            ///< Retrieve parameters.
    virtual void _setq(double*);           ///< Assign parameters.
    virtual void getqname(char**);         ///< Retrieve parameter names.

    virtual bool setChild(int index, Implicit* child);  ///< Assign elements.
    virtual Implicit* getChild(int index); ///< Retrieve elements.

    int numChildren() { return size(); }
    int maxChildren() { return -1; }    ///< Arbitrary # of children

    MAKE_NAME();
};

#endif

