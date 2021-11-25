pspg_tests_system_=pspg/tests/system/Test.cu

pspg_tests_system_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pspg_tests_system_))
pspg_tests_system_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspg_tests_system_:.cu=.o))

