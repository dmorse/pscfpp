#ifndef PRDC_TMPL_DECLARE_H
#define PRDC_TMPL_DECLARE_H

#ifdef PSCF_CPP
namespace Pscf {
   namespace Rpc {
      template <int D> class Types;
   }
}
#endif
#ifdef PSCF_CUDA
namespace Pscf {
   namespace Rpg {
      template <int D> class Types;
   }
}
#endif

// Output an input unchanged iff PSCF_CPP is defined
#ifdef PSCF_CPP
   #define IFDEF_PSCF_CPP(X) X
#else
   #define IFDEF_PSCF_CPP(X)
#endif

// Output an input unchanged iff PSCF_CUDA is defined
#ifdef PSCF_CUDA
   #define IFDEF_PSCF_CUDA(X) X
#else
   #define IFDEF_PSCF_CUDA(X)
#endif

// Declare explicit instantiations of all C++ specializations
#define PRDC_TMPL_DECLARE_CPP(CLASS_NAME) \
   extern template class CLASS_NAME <1, Rpc::Types<1> >; \
   extern template class CLASS_NAME <2, Rpc::Types<2> >; \
   extern template class CLASS_NAME <3, Rpc::Types<3> >; 

// Declare explicit instantiations of all CUDA specializations
#define PRDC_TMPL_DECLARE_CUDA(CLASS_NAME) \
   extern template class CLASS_NAME <1, Rpg::Types<1> >; \
   extern template class CLASS_NAME <2, Rpg::Types<2> >; \
   extern template class CLASS_NAME <3, Rpg::Types<3> >; 

// Declare explicit instantiations of all required specializations
#define PRDC_TMPL_DECLARE(CLASS_NAME) \
   IFDEF_PSCF_CPP( PRDC_TMPL_DECLARE_CPP(CLASS_NAME) ) \
   IFDEF_PSCF_CUDA( PRDC_TMPL_DECLARE_CUDA(CLASS_NAME) )

#endif // header guard
