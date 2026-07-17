#ifndef PRDC_TMPL_DEFINE_H
#define PRDC_TMPL_DEFINE_H

// Define explicit instantiations for all C++ specializations
#define PSCF_TMPL_DEFINE_CPP(CLASS_NAME) \
   template class CLASS_NAME <1, Rpc::Types<1> >; \
   template class CLASS_NAME <2, Rpc::Types<2> >; \
   template class CLASS_NAME <3, Rpc::Types<3> >; 

// Define explicit instantiations for all CUDA specializations
#define PSCF_TMPL_DEFINE_CUDA(CLASS_NAME) \
   template class CLASS_NAME <1, Rpg::Types<1> >; \
   template class CLASS_NAME <2, Rpg::Types<2> >; \
   template class CLASS_NAME <3, Rpg::Types<3> >; 

#endif
