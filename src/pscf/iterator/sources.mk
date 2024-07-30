pscf_iterator_= \
  pscf/iterator/AmbdInteraction.cpp \
  pscf/iterator/FieldGenerator.cpp \
  pscf/iterator/ImposedFieldsTmpl.cpp \
  pscf/iterator/NanException.cpp

pscf_iterator_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_iterator_:.cpp=.o))

