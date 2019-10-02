pspc_tests_field_=pspc/tests/field/Test.cc

pspc_tests_field_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pspc_tests_field_))
pspc_tests_field_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pspc_tests_field_:.cc=.o))

