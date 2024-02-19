r1d_domain_=\
  r1d/domain/GeometryMode.cpp \
  r1d/domain/Domain.cpp 

r1d_domain_OBJS=\
     $(addprefix $(BLD_DIR)/, $(r1d_domain_:.cpp=.o))

