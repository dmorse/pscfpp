pscf_tests_cuda_=pscf/tests/cuda/cudaTest.cu

pscf_tests_cuda_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_tests_cuda_:.cu=.o))

pscf_tests_cuda_DEPS=\
     $(addprefix $(BLD_DIR)/, $(pscf_tests_cuda_:.cu=.du))

