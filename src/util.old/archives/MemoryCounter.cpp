/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MemoryCounter.h"

namespace Util
{

   /*
   * Constructor.
   */
   MemoryCounter::MemoryCounter()
    : size_(0),
      version_(0)
   {}

   /*
   * Destructor.
   */
   MemoryCounter::~MemoryCounter()
   {}  

   /*
   * Return cursor to beginning.
   */
   void MemoryCounter::clear()
   {  size_ = 0; }

}
