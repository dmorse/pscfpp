#ifndef UTIL_BEGIN_H
#define UTIL_BEGIN_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComponent.h>
#include <util/param/Label.h>

#include <string>

namespace Util
{

   /**
   * Beginning line of a composite parameter block.
   *
   * \ingroup Param_Module
   */
   class Begin : public ParamComponent
   {

   public:

      /**
      * Constructor.
      */
      Begin(const char* label, bool isRequired = true);

      // Default destructor.

      /**
      * Read the opening line.
      *
      * \param in input stream
      */
      virtual void readParam(std::istream &in);

      /**
      * Write the opening line.
      *
      * \param out output stream
      */
      virtual void writeParam(std::ostream &out);

      /**
      * Is this the beginning line for a required element?
      */
      bool isRequired() const;

      /**
      * Is this an active element (has it been read from file)?
      */
      bool isActive() const;

      /**
      * Do-nothing implementation of virtual resetParam function.
      */
      virtual void resetParam();

   private:

      /// Classname label string (classname + "{")
      Label label_;

      /// Is this active (always true if isRequired).
      bool isActive_;

   };

   /*
   * Is this the beginning line for a required element?
   */
   inline bool Begin::isRequired() const
   {  return label_.isRequired(); }

   /*
   * Is this element active (has it been read from file)?
   */
   inline bool Begin::isActive() const
   {  return isActive_; }

} 
#endif
