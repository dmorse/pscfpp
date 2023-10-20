pscf_openmp_= \
  pscf/openmp/getNThread.cpp 

pscf_openmp_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_openmp_:.cpp=.o))

