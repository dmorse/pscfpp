#ifndef UTIL_INIT_STATIC_H
#define UTIL_INIT_STATIC_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace Util
{

   /**
   * Guarantee initialization of all static class members in Util namespace.
   */
   void initStatic();

}
#endif
