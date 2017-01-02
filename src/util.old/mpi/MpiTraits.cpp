/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MpiTraits.h"
#ifdef UTIL_MPI

namespace Util 
{

   const MPI::Datatype MpiTraitsNoType::type = MPI::CHAR; 
   const bool MpiTraitsNoType::hasType = false;

   const MPI::Datatype MpiTraits<char>::type = MPI::CHAR; 
   const bool MpiTraits<char>::hasType = true;

   const MPI::Datatype MpiTraits<unsigned char>::type = MPI::UNSIGNED_CHAR; 
   const bool MpiTraits<unsigned char>::hasType = true;

   const MPI::Datatype MpiTraits<short>::type = MPI::SHORT_INT; 
   const bool MpiTraits<short>::hasType = true;

   const MPI::Datatype MpiTraits<int>::type  = MPI::INT; 
   const bool MpiTraits<int>::hasType = true;

   const MPI::Datatype MpiTraits<long>::type = MPI::LONG; 
   const bool MpiTraits<long>::hasType = true;

   const MPI::Datatype MpiTraits<unsigned short>::type = MPI::UNSIGNED_SHORT; 
   const bool MpiTraits<unsigned short>::hasType = true;

   const MPI::Datatype MpiTraits<unsigned int>::type = MPI::UNSIGNED; 
   const bool MpiTraits<unsigned int>::hasType = true;

   const MPI::Datatype MpiTraits<unsigned long >::type = MPI::UNSIGNED_LONG; 
   const bool MpiTraits<unsigned long>::hasType = true;

   const MPI::Datatype MpiTraits<float>::type = MPI::FLOAT; 
   const bool MpiTraits<float>::hasType = true;

   const MPI::Datatype MpiTraits<double>::type = MPI::DOUBLE; 
   const bool MpiTraits<double>::hasType = true;

   const MPI::Datatype MpiTraits<long double>::type = MPI::LONG_DOUBLE; 
   const bool MpiTraits<long double>::hasType = true;

   const MPI::Datatype MpiTraits<bool>::type = MPI::BYTE; 
   const bool MpiTraits<bool>::hasType = false;

   #if 0
   const MPI::Datatype MpiTraits<std::complex<float> >::type = MPI::COMPLEX; 
   const bool MpiTraits<std::complex<float> >::hasType = true;

   const MPI::Datatype MpiTraits<std::complex<double> >::type = MPI::DOUBLE_COMPLEX; 
   const bool MpiTraits<std::complex<double> >::hasType = true;

   const MPI::Datatype MpiTraits<std::complex<long double> >::type = MPI::LONG_DOUBLE_COMPLEX; 
   const bool MpiTraits<std::complex<long double> >::hasType = true;
   #endif

   #if 0
   const MPI::Datatype MpiTraits<wchar_t>::type = MPI::WCHAR; 
   const bool MpiTraits<wchar_t>::hasType = true;
   #endif
}
#endif
