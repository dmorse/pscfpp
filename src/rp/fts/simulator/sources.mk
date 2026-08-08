
rp_fts_simulator_OBJS=
rp_fts_simulator_DEPS=

ifdef PSCF_CPP
  rp_fts_simulator_CPP= \
    rp/fts/simulator/Simulator.cpp \
    rp/fts/simulator/SimState.cpp  \
    rp/fts/simulator/SimulatorFactory.cpp
  rp_fts_simulator_OBJS+=\
       $(addprefix $(BLD_DIR)/, $(rp_fts_simulator_CPP:.cpp=.o))
  rp_fts_simulator_DEPS+=\
       $(addprefix $(BLD_DIR)/, $(rp_fts_simulator_CPP:.cpp=.d))
endif

ifdef PSCF_CUDA
  rp_fts_simulator_CU= \
    rp/fts/simulator/Simulator.cu \
    rp/fts/simulator/SimState.cu \
    rp/fts/simulator/SimulatorFactory.cu
  rp_fts_simulator_OBJS+=\
      $(addprefix $(BLD_DIR)/, $(rp_fts_simulator_CU:.cu=.ou))
  rp_fts_simulator_DEPS+=\
      $(addprefix $(BLD_DIR)/, $(rp_fts_simulator_CU:.cu=.du))
endif
