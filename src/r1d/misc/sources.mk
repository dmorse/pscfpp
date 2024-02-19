r1d_misc_=\
  r1d/misc/HomogeneousComparison.cpp \
  r1d/misc/FieldIo.cpp 

r1d_misc_OBJS=\
     $(addprefix $(BLD_DIR)/, $(r1d_misc_:.cpp=.o))

