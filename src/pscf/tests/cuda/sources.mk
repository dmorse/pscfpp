pscf_tests_cuda_=pscf/tests/cuda/Test.cu

pscf_tests_cuda_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pscf_tests_cuda_:.cu=.o))

