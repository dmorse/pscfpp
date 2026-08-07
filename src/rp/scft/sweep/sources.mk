
rp_scft_sweep_OBJS=
rp_scft_sweep_DEPS=

ifdef PSCF_CPP
  rp_scft_sweep_CPP= \
    rp/scft/sweep/BasisFieldState.cpp \
    rp/scft/sweep/Sweep.cpp \
    rp/scft/sweep/SweepParameter.cpp \
    rp/scft/sweep/LinearSweep.cpp \
    rp/scft/sweep/SweepFactory.cpp
  rp_scft_sweep_OBJS+= \
       $(addprefix $(BLD_DIR)/, $(rp_scft_sweep_CPP:.cpp=.o))
  rp_scft_sweep_DEPS+= \
       $(addprefix $(BLD_DIR)/, $(rp_scft_sweep_CPP:.cpp=.d))
endif

ifdef PSCF_CUDA
  rp_scft_sweep_CU= \
    rp/scft/sweep/BasisFieldState.cu \
    rp/scft/sweep/Sweep.cu \
    rp/scft/sweep/SweepParameter.cu \
    rp/scft/sweep/LinearSweep.cu \
    rp/scft/sweep/SweepFactory.cu
  rp_scft_sweep_OBJS+= \
       $(addprefix $(BLD_DIR)/, $(rp_scft_sweep_CU:.cu=.ou))
  rp_scft_sweep_DEPS+= \
       $(addprefix $(BLD_DIR)/, $(rp_scft_sweep_CU:.cu=.du))
endif

