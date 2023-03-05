cyln_tests_field_=cyln/tests/field/Test.cc

cyln_tests_field_SRCS=\
     $(addprefix $(SRC_DIR)/, $(cyln_tests_field_))
cyln_tests_field_OBJS=\
     $(addprefix $(BLD_DIR)/, $(cyln_tests_field_:.cc=.o))

