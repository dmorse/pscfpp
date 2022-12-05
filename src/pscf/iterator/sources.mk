pscf_iterator_= \
  pscf/iterator/AmbdResidual.cpp 

pscf_iterator_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pscf_iterator_))
pscf_iterator_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_iterator_:.cpp=.o))

