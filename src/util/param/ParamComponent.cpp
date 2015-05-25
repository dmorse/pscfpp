/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ParamComponent.h"

namespace Util
{

   bool ParamComponent::echo_ = false;

   /*
   * Constructor.
   *
   * Note: Default constructor for _indent string creates an empty string.
   */
   ParamComponent::ParamComponent() 
    : MpiFileIo(),
      indent_()
   {}

   /*
   * Copy constructor.
   */
   ParamComponent::ParamComponent(const ParamComponent& other) 
    : MpiFileIo(other),
      indent_(other.indent_)
   {}

   /*
   * Destructor.
   */
   ParamComponent::~ParamComponent() 
   {} 

   /*
   * Return indent string for this object (string of spaces).
   */
   std::string ParamComponent::indent() const
   {  return indent_; }

   /*
   * Set indent level to be one higher than that of parent.
   */
   void ParamComponent::setIndent(const ParamComponent& parent, bool next)
   { 
      indent_ = parent.indent();
      if (next) {
         std::string space("  ");
         indent_ += space; 
      }
   }

   // Static functions

   /*
   * Set echo to true (default) or false.
   */
   void ParamComponent::setEcho(bool echo)
   {  echo_ = echo; }

   /*
   * Get value of echo. true = echoing is on, false = off.
   */
   bool ParamComponent::echo()
   {  return echo_; }

   /*
   * This static method exists to guarantee initialization of a static 
   * constant echo_ that is defined in the same file.  Call it somewhere 
   * in the program to guarantee that the contents of this file will 
   * be linked, rather than optimized away. It may only be called once.
   */
   void ParamComponent::initStatic()
   {  
      // This function can only be called once.
      static int nCall = 0;
      if (nCall == 0) {
         echo_ = false; 
      }
      ++nCall;
   }

} 
