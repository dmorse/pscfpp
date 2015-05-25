#ifndef UTIL_BLANK_H
#define UTIL_BLANK_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComponent.h>

namespace Util
{

   /**
   * An empty line within a parameter file. 
   *
   * A Param represents an empty line within a file format that is
   * represented as a ParamComposite.
   *
   * \ingroup Param_Module
   */
   class Blank : public ParamComponent
   {

   public:

      /// Constructor.
      Blank();

      /// Virtual Destructor
      virtual ~Blank();

      /**
      * Read a blank line
      *
      * \param in input stream
      */
      virtual void readParam(std::istream &in);

      /**
      * Write a blank line
      *
      * \param out output stream
      */
      virtual void writeParam(std::ostream &out);

   };

} 
#endif
