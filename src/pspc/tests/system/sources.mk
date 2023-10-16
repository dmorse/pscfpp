pspc_tests_system_=pspc/tests/system/Test.cpp

pspc_tests_system_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspc_tests_system_:.cpp=.o))

