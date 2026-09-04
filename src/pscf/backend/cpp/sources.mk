pscf_backend_cpp_CPP= \
  pscf/backend/cpp/complex.cpp \
  pscf/backend/cpp/VecOp.cpp \
  pscf/backend/cpp/VecOpCx.cpp \
  pscf/backend/cpp/Reduce.cpp \
  pscf/backend/cpp/ReduceCx.cpp \
  pscf/backend/cpp/CpuVecRandom.cpp 

pscf_backend_cpp_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_backend_cpp_CPP:.cpp=.o))

pscf_backend_cpp_DEPS=\
     $(addprefix $(BLD_DIR)/, $(pscf_backend_cpp_CPP:.cpp=.d))
