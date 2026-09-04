#ifndef PRDC_BACKEND_ID_H
#define PRDC_BACKEND_ID_H

/*
* PSCF - Polymer Self-Consistent Field 
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace Pscf {
namespace Prdc {

   /**
   * Enum that identifies a computational "backend".
   *
   * This enum is used in template code that can use any of several 
   * different computational "backends" to perform basic mathematical 
   * operations that some backends can perform in parallel.
   * 
   * Enum values:
   *
   *    - Cpp  : serial Cpu backend, programmed in standard C++
   *    - Cuda : Gpu backend, programmed in C++ Cuda
   *    - None : default null value
   */
   enum class BackendId {
      Cpp, Cuda, None
   };


} // namespace Prdc
} // namespace Pscf
#endif
