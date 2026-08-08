
rp_fts_brownian_OBJS=
rp_fts_brownian_DEPS=

ifdef PSCF_CPP
  rp_fts_brownian_CPP= \
    rp/fts/brownian/BdStep.cpp \
    rp/fts/brownian/BdSimulator.cpp \
    rp/fts/brownian/ExplicitBdStep.cpp \
    rp/fts/brownian/PredCorrBdStep.cpp \
    rp/fts/brownian/LMBdStep.cpp \
    rp/fts/brownian/BdStepFactory.cpp
  rp_fts_brownian_OBJS+=\
       $(addprefix $(BLD_DIR)/, $(rp_fts_brownian_CPP:.cpp=.o))
  rp_fts_brownian_DEPS+=\
       $(addprefix $(BLD_DIR)/, $(rp_fts_brownian_CPP:.cpp=.d))
endif

ifdef PSCF_CUDA
  rp_fts_brownian_CU= \
    rp/fts/brownian/BdStep.cu \
    rp/fts/brownian/BdSimulator.cu \
    rp/fts/brownian/ExplicitBdStep.cu \
    rp/fts/brownian/PredCorrBdStep.cu \
    rp/fts/brownian/LMBdStep.cu \
    rp/fts/brownian/BdStepFactory.cu
  rp_fts_brownian_OBJS+=\
      $(addprefix $(BLD_DIR)/, $(rp_fts_brownian_CU:.cu=.ou))
  rp_fts_brownian_DEPS+=\
      $(addprefix $(BLD_DIR)/, $(rp_fts_brownian_CU:.cu=.du))
endif
