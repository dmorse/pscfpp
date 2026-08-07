rp_solvers_OBJS=
rp_solvers_DEPS=

ifdef PSCF_CPP
  rp_solvers_CPP= \
    rp/solvers/Propagator_c.cpp \
    rp/solvers/Block_c.cpp \
    rp/solvers/Polymer_c.cpp \
    rp/solvers/Solvent_c.cpp \
    rp/solvers/Mixture_c.cpp \
    rp/solvers/MixtureModifier_c.cpp
  rp_solvers_OBJS+=\
     $(addprefix $(BLD_DIR)/, $(rp_solvers_CPP:.cpp=.o))
  rp_solvers_DEPS+=\
     $(addprefix $(BLD_DIR)/, $(rp_solvers_CPP:.cpp=.d))
endif

ifdef PSCF_CUDA
  rp_solvers_CU= \
    rp/solvers/Propagator_u.cu \
    rp/solvers/Block_u.cu \
    rp/solvers/Polymer_u.cu \
    rp/solvers/Solvent_u.cu \
    rp/solvers/Mixture_u.cu \
    rp/solvers/MixtureModifier_u.cu
  rp_solvers_OBJS+=\
       $(addprefix $(BLD_DIR)/, $(rp_solvers_CU:.cu=.o))
  rp_solvers_DEPS+=\
     $(addprefix $(BLD_DIR)/, $(rp_solvers_CU:.cu=.du))
endif

