
rp_fts_perturbation_OBJS=
rp_fts_perturbation_DEPS=

ifdef PSCF_CPP
  rp_fts_perturbation_CPP= \
    rp/fts/perturbation/Perturbation.cpp \
    rp/fts/perturbation/EinsteinCrystalPerturbation.cpp \
    rp/fts/perturbation/PerturbationFactory.cpp
  rp_fts_perturbation_OBJS+=\
       $(addprefix $(BLD_DIR)/, $(rp_fts_perturbation_CPP:.cpp=.o))
  rp_fts_perturbation_DEPS+=\
       $(addprefix $(BLD_DIR)/, $(rp_fts_perturbation_CPP:.cpp=.d))
endif

ifdef PSCF_CUDA
  rp_fts_perturbation_CU= \
    rp/fts/perturbation/Perturbation.cu \
    rp/fts/perturbation/EinsteinCrystalPerturbation.cu \
    rp/fts/perturbation/PerturbationFactory.cu
  rp_fts_perturbation_OBJS+=\
      $(addprefix $(BLD_DIR)/, $(rp_fts_perturbation_CU:.cu=.ou))
  rp_fts_perturbation_DEPS+=\
      $(addprefix $(BLD_DIR)/, $(rp_fts_perturbation_CU:.cu=.du))
endif
