fd1d_domain_=\
  fd1d/domain/GeometryMode.cpp \
  fd1d/domain/Domain.cpp 

fd1d_domain_OBJS=\
     $(addprefix $(BLD_DIR)/, $(fd1d_domain_:.cpp=.o))

