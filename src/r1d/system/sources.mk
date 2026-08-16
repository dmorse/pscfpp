r1d_system_CPP=\
  r1d/system/System.cpp \
  r1d/system/SystemAccess.cpp \
  r1d/system/HomogeneousComparison.cpp 

r1d_system_OBJS=\
     $(addprefix $(BLD_DIR)/, $(r1d_system_CPP:.cpp=.o))

r1d_system_DEPS=\
     $(addprefix $(BLD_DIR)/, $(r1d_system_CPP:.cpp=.d))

