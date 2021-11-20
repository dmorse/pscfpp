pspc_tests_sweep_=pspc/tests/sweep/Test.cc

pspc_tests_sweep_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pspc_tests_sweep_))
pspc_tests_sweep_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspc_tests_sweep_:.cc=.o))

