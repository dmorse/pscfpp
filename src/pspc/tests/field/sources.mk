pspc_tests_field_=pspc/tests/field/Test.cpp

pspc_tests_field_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspc_tests_field_:.cpp=.o))

