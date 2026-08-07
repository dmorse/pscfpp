ifdef PSCF_CUDA
  pscf_tests_=pscf/tests/cudaTest.cu
  pscf_tests_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_tests_:.cu=.ou))
else
  pscf_tests_=pscf/tests/cpuTest.cpp
  pscf_tests_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_tests_:.cpp=.o))
endif

