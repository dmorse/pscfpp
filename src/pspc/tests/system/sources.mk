pspc_tests_system_=pspc/tests/system/Test.cc

pspc_tests_system_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspc_tests_system_:.cc=.o))

