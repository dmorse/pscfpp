
rp_fts_montecarlo_OBJS=
rp_fts_montecarlo_DEPS=

ifdef PSCF_CPP
  rp_fts_montecarlo_CPP= \
    rp/fts/montecarlo/McMove.cpp \
    rp/fts/montecarlo/McMoveManager.cpp \
    rp/fts/montecarlo/McSimulator.cpp \
    rp/fts/montecarlo/RealMove.cpp \
    rp/fts/montecarlo/ForceBiasMove.cpp \
    rp/fts/montecarlo/BdMove.cpp \
    rp/fts/montecarlo/ShiftMove.cpp \
    rp/fts/montecarlo/McMoveFactory.cpp
  rp_fts_montecarlo_OBJS+=\
       $(addprefix $(BLD_DIR)/, $(rp_fts_montecarlo_CPP:.cpp=.o))
  rp_fts_montecarlo_DEPS+=\
       $(addprefix $(BLD_DIR)/, $(rp_fts_montecarlo_CPP:.cpp=.d))
endif

ifdef PSCF_CUDA
  rp_fts_montecarlo_CU= \
    rp/fts/montecarlo/McMove.cu \
    rp/fts/montecarlo/McMoveManager.cu \
    rp/fts/montecarlo/McSimulator.cu \
    rp/fts/montecarlo/RealMove.cu \
    rp/fts/montecarlo/ForceBiasMove.cu \
    rp/fts/montecarlo/BdMove.cu \
    rp/fts/montecarlo/ShiftMove.cu \
    rp/fts/montecarlo/McMoveFactory.cu
  rp_fts_montecarlo_OBJS+=\
      $(addprefix $(BLD_DIR)/, $(rp_fts_montecarlo_CU:.cu=.ou))
  rp_fts_montecarlo_DEPS+=\
      $(addprefix $(BLD_DIR)/, $(rp_fts_montecarlo_CU:.cu=.du))
endif
