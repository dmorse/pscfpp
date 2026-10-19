#ifndef PSCF_TMPL_DEFINE_H
#define PSCF_TMPL_DEFINE_H

// Define explicit instantiations for all C++ specializations
#define PSCF_TMPL_DEFINE_CPP(CLASS_NAME) \
   template class CLASS_NAME <1, CPT >; \
   template class CLASS_NAME <2, CPT >; \
   template class CLASS_NAME <3, CPT >; 

// Define explicit instantiations for all CUDA specializations
#define PSCF_TMPL_DEFINE_CUDA(CLASS_NAME) \
   template class CLASS_NAME <1, CUT >; \
   template class CLASS_NAME <2, CUT >; \
   template class CLASS_NAME <3, CUT >; 

#endif
