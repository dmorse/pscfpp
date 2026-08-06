pscf_cpu_CPP= \
  pscf/cpu/complex.cpp \
  pscf/cpu/VecOp.cpp \
  pscf/cpu/VecOpCx.cpp \
  pscf/cpu/Reduce.cpp \
  pscf/cpu/ReduceCx.cpp \
  pscf/cpu/CpuVecRandom.cpp 

pscf_cpu_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_cpu_CPP:.cpp=.o))

