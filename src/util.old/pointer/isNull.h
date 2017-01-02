#ifndef UTIL_IS_NULL_H
#define UTIL_IS_NULL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace Util
{

   /** 
   * Return true iff a built-in pointer is null.
   */
   template <typename T>
   inline bool isNull(T* ptr)
   {  return (ptr == 0); }

}
#endif 
