pssp_tests_system_=pssp/tests/system/Test.cc

pssp_tests_system_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pssp_tests_system_))
pssp_tests_system_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pssp_tests_system_:.cc=.o))

