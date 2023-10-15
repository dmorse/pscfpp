prdc_tests_=prdc/tests/Test.cu

prdc_tests_SRCS=\
     $(addprefix $(SRC_DIR)/, $(prdc_tests_))
prdc_tests_OBJS=\
     $(addprefix $(BLD_DIR)/, $(prdc_tests_:.cu=.o))

