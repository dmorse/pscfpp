r1_system_CPP=\
  r1/system/System.cpp \
  r1/system/SystemAccess.cpp \
  r1/system/HomogeneousComparison.cpp 

r1_system_OBJS=\
     $(addprefix $(BLD_DIR)/, $(r1_system_CPP:.cpp=.o))

r1_system_DEPS=\
     $(addprefix $(BLD_DIR)/, $(r1_system_CPP:.cpp=.d))

