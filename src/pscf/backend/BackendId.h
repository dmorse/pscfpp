#ifndef PSCF_BACKEND_ID_H
#define PSCF_BACKEND_ID_H

/*
* PSCF - Polymer Self-Consistent Field 
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

namespace Pscf {

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
   */
   enum class BackendId {
      Cpp, Cuda
   };


} // namespace Pscf
#endif
