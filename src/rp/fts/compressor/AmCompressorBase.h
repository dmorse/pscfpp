#ifndef RP_AM_COMPRESSOR_BASE_H
#define RP_AM_COMPRESSOR_BASE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/iterator/AmIteratorTmpl.h>  // base class template
#include <rp/fts/compressor/Compressor.h>  // template argument
#include <pscf/backend/BackendArray.h>      // template argument

// Explicit instantiation declarations
namespace Pscf {
   #ifdef PSCF_CPP
   class CPT;
   extern template 
   class AmIteratorTmpl< Rp::Compressor<1,CPT>, BackendArray<double,CPT> >;
   extern template 
   class AmIteratorTmpl< Rp::Compressor<2,CPT>, BackendArray<double,CPT> >; 
   extern template 
   class AmIteratorTmpl< Rp::Compressor<3,CPT>, BackendArray<double,CPT> >;
   #endif
   #ifdef PSCF_CUDA
   class CUT;
   extern template 
   class AmIteratorTmpl< Rp::Compressor<1,CUT>, BackendArray<double,CUT> >;
   extern template 
   class AmIteratorTmpl< Rp::Compressor<2,CUT>, BackendArray<double,CUT> >; 
   extern template 
   class AmIteratorTmpl< Rp::Compressor<3,CUT>, BackendArray<double,CUT> >;
   #endif
}

#endif
