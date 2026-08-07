ifdef PSCF_CUDA
  prdc_tests_=prdc/tests/cudaTest.cu
  prdc_tests_OBJS=\
     $(addprefix $(BLD_DIR)/, $(prdc_tests_:.cu=.ou))
  prdc_tests_DEPS=\
     $(addprefix $(BLD_DIR)/, $(prdc_tests_:.cu=.du))
else
  prdc_tests_=prdc/tests/cpuTest.cpp
  prdc_tests_OBJS=\
     $(addprefix $(BLD_DIR)/, $(prdc_tests_:.cpp=.o))
  prdc_tests_DEPS=\
     $(addprefix $(BLD_DIR)/, $(prdc_tests_:.cpp=.d))
endif

