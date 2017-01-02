#ifndef UTIL_OFFSET_H
#define UTIL_OFFSET_H

#include <cstddef>

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace Util
{

   /**
   * Template for calculating offsets of data members.
   *
   * Types: 
   *  D - derived class
   *  B - base class 
   *  M - member type
   * 
   * \ingroup Misc_Module
   */
   template <typename D, typename B, typename M>
   ptrdiff_t memberOffset(D& object, M  B::* memPtr)
   {  
      return reinterpret_cast<char*>(&(object.*memPtr)) - reinterpret_cast<char*>(&object); 
   }
   
   /**
   * Template for calculating offsets of base class subobjects.
   *
   * Types: 
   *  D - derived class
   *  B - base class 
   */
   template <typename D, typename B>
   ptrdiff_t baseOffset(D& object)
   {  
      D* d = &object; 
      B* b = (B*) d;
      return reinterpret_cast<char*>(d) - reinterpret_cast<char*>(b); 
   }
   


}
#endif
