#ifndef PSCF_TMPL_DECLARE_H
#define PSCF_TMPL_DECLARE_H

// Include required backend type classes
#ifdef PSCF_CPP
#include <pscf/backend/cpp/CPT.h>
#endif
#ifdef PSCF_CUDA
#include <pscf/backend/cuda/CUT.h>
#endif

// Return an input unchanged ifdef PSCF_CPP, or nothing otherwise
#ifdef PSCF_CPP
   #define IFDEF_PSCF_CPP(X) X
#else
   #define IFDEF_PSCF_CPP(X)
#endif

// Return an input unchanged ifdef PSCF_CUDA, or nothing otherwise
#ifdef PSCF_CUDA
   #define IFDEF_PSCF_CUDA(X) X
#else
   #define IFDEF_PSCF_CUDA(X)
#endif

// Declare explicit instantiations of all C++ specializations
#define PSCF_TMPL_DECLARE_CPP(CLASS_NAME) \
   extern template class CLASS_NAME <1,CPT>; \
   extern template class CLASS_NAME <2,CPT>; \
   extern template class CLASS_NAME <3,CPT>; 

// Declare explicit instantiations of all CUDA specializations
#define PSCF_TMPL_DECLARE_CUDA(CLASS_NAME) \
   extern template class CLASS_NAME <1,CUT>; \
   extern template class CLASS_NAME <2,CUT>; \
   extern template class CLASS_NAME <3,CUT>; 

// Declare explicit instantiations of all required specializations
#define PSCF_TMPL_DECLARE(CLASS_NAME) \
   IFDEF_PSCF_CPP( PSCF_TMPL_DECLARE_CPP(CLASS_NAME) ) \
   IFDEF_PSCF_CUDA( PSCF_TMPL_DECLARE_CUDA(CLASS_NAME) )

#endif // header guard
