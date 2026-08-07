
rp_scft_iterator_OBJS=
rp_scft_iterator_DEPS=

ifdef PSCF_CPP
  rp_scft_iterator_CPP= \
    rp/scft/iterator/Iterator.cpp \
    rp/scft/iterator/AmIteratorBasis.cpp \
    rp/scft/iterator/AmIteratorGrid.cpp \
    rp/scft/iterator/IteratorFactory.cpp
  rp_scft_iterator_OBJS+=\
       $(addprefix $(BLD_DIR)/, $(rp_scft_iterator_CPP:.cpp=.o))
  rp_scft_iterator_DEPS+=\
       $(addprefix $(BLD_DIR)/, $(rp_scft_iterator_CPP:.cpp=.d))
endif

ifdef PSCF_CUDA
  rp_scft_iterator_CU= \
    rp/scft/iterator/Iterator.cu \
    rp/scft/iterator/AmIteratorBasis.cu \
    rp/scft/iterator/AmIteratorGrid.cu \
    rp/scft/iterator/IteratorFactory.cu
  rp_scft_iterator_OBJS+=\
      $(addprefix $(BLD_DIR)/, $(rp_scft_iterator_CU:.cu=.ou))
  rp_scft_iterator_DEPS+=\
      $(addprefix $(BLD_DIR)/, $(rp_scft_iterator_CU:.cu=.du))
endif
