cpc_tests_field_=cpc/tests/field/Test.cpp

cpc_tests_field_OBJS=\
     $(addprefix $(BLD_DIR)/, $(cpc_tests_field_:.cpp=.o))

