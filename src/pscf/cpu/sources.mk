pscf_cpu_= \
  pscf/cpu/complex.cpp \
  pscf/cpu/VecOp.cpp \
  pscf/cpu/VecOpCx.cpp \
  pscf/cpu/Reduce.cpp \
  pscf/cpu/ReduceCx.cpp \
  pscf/cpu/CpuVecRandom.cpp \
  pscf/cpu/CppTp.cpp

pscf_cpu_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_cpu_:.cpp=.o))

