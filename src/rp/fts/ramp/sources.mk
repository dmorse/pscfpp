
rp_fts_ramp_OBJS=
rp_fts_ramp_DEPS=

ifdef PSCF_CPP
  rp_fts_ramp_CPP= \
    rp/fts/ramp/Ramp.cpp \
    rp/fts/ramp/RampParameter.cpp \
    rp/fts/ramp/LinearRamp.cpp \
    rp/fts/ramp/RampFactory.cpp
  rp_fts_ramp_OBJS+=\
       $(addprefix $(BLD_DIR)/, $(rp_fts_ramp_CPP:.cpp=.o))
  rp_fts_ramp_DEPS+=\
       $(addprefix $(BLD_DIR)/, $(rp_fts_ramp_CPP:.cpp=.d))
endif

ifdef PSCF_CUDA
  rp_fts_ramp_CU= \
    rp/fts/ramp/Ramp.cu \
    rp/fts/ramp/RampParameter.cu \
    rp/fts/ramp/LinearRamp.cu \
    rp/fts/ramp/RampFactory.cu
  rp_fts_ramp_OBJS+=\
      $(addprefix $(BLD_DIR)/, $(rp_fts_ramp_CU:.cu=.ou))
  rp_fts_ramp_DEPS+=\
      $(addprefix $(BLD_DIR)/, $(rp_fts_ramp_CU:.cu=.du))
endif
