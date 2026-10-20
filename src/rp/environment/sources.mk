rp_environment_OBJS_=
rp_environment_DEPS_=

rp_environment_CPP= \
  rp/environment/Environment.cpp \
  rp/environment/FieldGenerator.cpp \
  rp/environment/MixAndMatchEnv.cpp \
  rp/environment/FilmFieldGenMaskBase.cpp \
  rp/environment/FilmFieldGenExtBase.cpp

ifdef PSCF_CPP
  rp_environment_CPP+= \
     rp/environment/FilmFieldGenExt.cpp \
     rp/environment/FilmFieldGenMask.cpp \
     rp/environment/FilmEnvironment.cpp \
     rp/environment/EnvironmentFactory.cpp
endif
rp_environment_OBJS+= \
       $(addprefix $(BLD_DIR)/, $(rp_environment_CPP:.cpp=.o))
rp_environment_DEPS+= \
       $(addprefix $(BLD_DIR)/, $(rp_environment_CPP:.cpp=.d))

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
