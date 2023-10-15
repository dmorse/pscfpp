prdc_tests_cuda_=prdc/tests/cuda/Test.ccu

prdc_tests_cuda_OBJS=\
     $(addprefix $(BLD_DIR)/, $(prdc_tests_cuda_:.ccu=.o))

