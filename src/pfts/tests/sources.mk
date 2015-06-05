pfts_tests_=pfts/tests/Test.cc

pfts_tests_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pfts_tests_))
pfts_tests_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pfts_tests_:.cc=.o))

