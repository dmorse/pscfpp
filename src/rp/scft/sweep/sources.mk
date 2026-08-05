
rp_scft_sweep_OBJS=
rp_scft_sweep_DEPS=

ifdef PSCF_CPP
  rp_scft_sweep_CPP= \
    rp/scft/sweep/BasisFieldState_c.cpp \
    rp/scft/sweep/Sweep_c.cpp \
    rp/scft/sweep/SweepParameter_c.cpp \
    rp/scft/sweep/LinearSweep_c.cpp \
    rp/scft/sweep/SweepFactory_c.cpp
  rp_scft_sweep_OBJS+= \
       $(addprefix $(BLD_DIR)/, $(rp_scft_sweep_CPP:.cpp=.o))
  rp_scft_sweep_DEPS+= \
       $(addprefix $(BLD_DIR)/, $(rp_scft_sweep_CPP:.cpp=.d))
endif

ifdef PSCF_CUDA
  rp_scft_sweep_CU= \
    rp/scft/sweep/BasisFieldState_u.cu \
    rp/scft/sweep/Sweep_u.cu \
    rp/scft/sweep/SweepParameter_u.cu \
    rp/scft/sweep/LinearSweep_u.cu \
    rp/scft/sweep/SweepFactory_u.cu
  rp_scft_sweep_OBJS+= \
       $(addprefix $(BLD_DIR)/, $(rp_scft_sweep_CU:.cu=.o))
  rp_scft_sweep_DEPS+= \
       $(addprefix $(BLD_DIR)/, $(rp_scft_sweep_CU:.cu=.d))
endif

