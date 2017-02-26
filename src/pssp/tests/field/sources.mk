pssp_tests_field_=pssp/tests/field/Test.cc

pssp_tests_field_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pssp_tests_field_))
pssp_tests_field_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pssp_tests_field_:.cc=.o))

