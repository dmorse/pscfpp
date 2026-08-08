rp_environment_OBJS_=
rp_environment_DEPS_=

ifdef PSCF_CPP
  rp_environment_CPP= \
     rp/environment/FilmFieldGenExt.cpp \
     rp/environment/FilmFieldGenMask.cpp \
     rp/environment/FilmEnvironment.cpp \
     rp/environment/EnvironmentFactory.cpp
  rp_environment_OBJS+= \
       $(addprefix $(BLD_DIR)/, $(rp_environment_CPP:.cpp=.o))
  rp_environment_DEPS+= \
       $(addprefix $(BLD_DIR)/, $(rp_environment_CPP:.cpp=.d))
endif

ifdef PSCF_CUDA
  rp_environment_CU= \
     rp/environment/FilmFieldGenExt.cu \
     rp/environment/FilmFieldGenMask.cu \
     rp/environment/FilmEnvironment.cu \
     rp/environment/EnvironmentFactory.cu
  rp_environment_OBJS+= \
     $(addprefix $(BLD_DIR)/, $(rp_environment_CU:.cu=.ou))
  rp_environment_DEPS+= \
     $(addprefix $(BLD_DIR)/, $(rp_environment_CU:.cu=.du))
endif
