prdc_tests_cuda_=prdc/tests/cuda/Test.cu

prdc_tests_cuda_OBJS=\
     $(addprefix $(BLD_DIR)/, $(prdc_tests_cuda_:.cu=.o))

prdc_tests_cuda_DEPS=\
     $(addprefix $(BLD_DIR)/, $(prdc_tests_cuda_:.cu=.d))

