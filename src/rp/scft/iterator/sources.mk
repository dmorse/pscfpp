
rp_scft_iterator_OBJS=

ifdef PSCF_CPP
rp_scft_iterator_CPP= \
  rp/scft/iterator/Iterator.cpp \
  rp/scft/iterator/AmIteratorBasis.cpp \
  rp/scft/iterator/AmIteratorGrid.cpp \
  rp/scft/iterator/IteratorFactory.cpp
rp_scft_iterator_OBJS+=\
     $(addprefix $(BLD_DIR)/, $(rp_scft_iterator_CPP:.cpp=.o))
endif

ifdef PSCF_CUDA
rp_scft_iterator_CUDA= \
  rp/scft/iterator/Iterator.cu \
  rp/scft/iterator/AmIteratorBasis.cu \
  rp/scft/iterator/AmIteratorGrid.cu \
  rp/scft/iterator/IteratorFactory.cu
rp_scft_iterator_OBJS+=\
     $(addprefix $(BLD_DIR)/, $(rp_scft_iterator_CUDA:.cu=.o))
endif
