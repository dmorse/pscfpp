rp_solvers_OBJS=
rp_solvers_DEPS=

ifdef PSCF_CPP
  rp_solvers_CPP= \
    rp/solvers/Propagator.cpp \
    rp/solvers/Block.cpp \
    rp/solvers/Polymer.cpp \
    rp/solvers/Solvent.cpp \
    rp/solvers/Mixture.cpp \
    rp/solvers/MixtureModifier.cpp
  rp_solvers_OBJS+=\
     $(addprefix $(BLD_DIR)/, $(rp_solvers_CPP:.cpp=.o))
  rp_solvers_DEPS+=\
     $(addprefix $(BLD_DIR)/, $(rp_solvers_CPP:.cpp=.d))
endif

ifdef PSCF_CUDA
  rp_solvers_CU= \
    rp/solvers/Propagator.cu \
    rp/solvers/Block.cu \
    rp/solvers/Polymer.cu \
    rp/solvers/Solvent.cu \
    rp/solvers/Mixture.cu \
    rp/solvers/MixtureModifier.cu
  rp_solvers_OBJS+=\
       $(addprefix $(BLD_DIR)/, $(rp_solvers_CU:.cu=.ou))
  rp_solvers_DEPS+=\
     $(addprefix $(BLD_DIR)/, $(rp_solvers_CU:.cu=.du))
endif

