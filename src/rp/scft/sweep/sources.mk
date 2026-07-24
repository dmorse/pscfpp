
rpc_scft_sweep_OBJS=

ifdef PSCF_CPP
rpc_scft_sweep_CPP= \
  rpc/scft/sweep/BasisFieldState_c.cpp \
  rpc/scft/sweep/Sweep_c.cpp \
  rpc/scft/sweep/SweepParameter_c.cpp \
  rpc/scft/sweep/LinearSweep_c.cpp \
  rpc/scft/sweep/SweepFactory_c.cpp
rpc_scft_sweep_OBJS+= \
     $(addprefix $(BLD_DIR)/, $(rpc_scft_sweep_CPP:.cpp=.o))
endif

ifdef PSCF_CUDA
rpc_scft_sweep_CUDA= \
  rpc/scft/sweep/BasisFieldState_u.cu \
  rpc/scft/sweep/Sweep_u.cu \
  rpc/scft/sweep/SweepParameter_u.cu \
  rpc/scft/sweep/LinearSweep_u.cu \
  rpc/scft/sweep/SweepFactory_u.cu
rpc_scft_sweep_OBJS+= \
     $(addprefix $(BLD_DIR)/, $(rpc_scft_sweep_CUDA:.cu=.o))
endif

