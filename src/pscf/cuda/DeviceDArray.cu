/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "DeviceDArray.tpp"

namespace Pscf {

   template class DeviceDArray<cudaReal>;
   template class DeviceDArray<cudaComplex>;

}