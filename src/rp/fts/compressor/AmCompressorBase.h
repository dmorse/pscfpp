#ifndef RP_AM_COMPRESSOR_BASE_H
#define RP_AM_COMPRESSOR_BASE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/iterator/AmIteratorTmpl.h>  // base class template
#include <rp/fts/compressor/Compressor.h>  // template argument
#include <pscf/backend/DeviceArray.h>      // template argument

// Explicit instantiation declarations
namespace Pscf {
   #ifdef PSCF_CPP
   class CPT;
   extern template 
   class AmIteratorTmpl< Rp::Compressor<1,CPT>, DeviceArray<double,CPT> >;
   extern template 
   class AmIteratorTmpl< Rp::Compressor<2,CPT>, DeviceArray<double,CPT> >; 
   extern template 
   class AmIteratorTmpl< Rp::Compressor<3,CPT>, DeviceArray<double,CPT> >;
   #endif
   #ifdef PSCF_CUDA
   class CUT;
   extern template 
   class AmIteratorTmpl< Rp::Compressor<1,CUT>, DeviceArray<double,CUT> >;
   extern template 
   class AmIteratorTmpl< Rp::Compressor<2,CUT>, DeviceArray<double,CUT> >; 
   extern template 
   class AmIteratorTmpl< Rp::Compressor<3,CUT>, DeviceArray<double,CUT> >;
   #endif
}

#endif
