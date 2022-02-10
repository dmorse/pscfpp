pspg_tests_cuda_=pspg/tests/cuda/Test.ccu

pspg_tests_cuda_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pspg_tests_cuda_))
pspg_tests_cuda_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspg_tests_cuda_:.ccu=.o))

