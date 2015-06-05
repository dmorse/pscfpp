pfts_tests_chem_=pfts/tests/chem/Test.cc

pfts_tests_chem_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pfts_tests_chem_))
pfts_tests_chem_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pfts_tests_chem_:.cc=.o))

