
rp_scft_iterator_OBJS=
rp_scft_iterator_DEPS=

ifdef PSCF_CPP
  rp_scft_iterator_CPP= \
    rp/scft/iterator/Iterator_c.cpp \
    rp/scft/iterator/AmIteratorBasis_c.cpp \
    rp/scft/iterator/AmIteratorGrid_c.cpp \
    rp/scft/iterator/IteratorFactory_c.cpp
  rp_scft_iterator_OBJS+=\
       $(addprefix $(BLD_DIR)/, $(rp_scft_iterator_CPP:.cpp=.o))
  rp_scft_iterator_DEPS+=\
       $(addprefix $(BLD_DIR)/, $(rp_scft_iterator_CPP:.cpp=.d))
endif

ifdef PSCF_CUDA
  rp_scft_iterator_CU= \
    rp/scft/iterator/Iterator_u.cu \
    rp/scft/iterator/AmIteratorBasis.cu \
    rp/scft/iterator/AmIteratorGrid.cu \
    rp/scft/iterator/IteratorFactory.cu
  rp_scft_iterator_OBJS+=\
      $(addprefix $(BLD_DIR)/, $(rp_scft_iterator_CU:.cu=.o))
  rp_scft_iterator_DEPS+=\
      $(addprefix $(BLD_DIR)/, $(rp_scft_iterator_CU:.cu=.d))
endif
