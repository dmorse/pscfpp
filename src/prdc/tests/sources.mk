ifdef PSCF_CUDA
  prdc_tests_=prdc/tests/cudaTest.cu
  prdc_tests_OBJS=\
     $(addprefix $(BLD_DIR)/, $(prdc_tests_:.cu=.o))
else
  prdc_tests_=prdc/tests/cpuTest.cpp
  prdc_tests_OBJS=\
     $(addprefix $(BLD_DIR)/, $(prdc_tests_:.cpp=.o))
endif

