chem_tests_=chem/tests/Test.cc

chem_tests_SRCS=\
     $(addprefix $(SRC_DIR)/, $(chem_tests_))
chem_tests_OBJS=\
     $(addprefix $(BLD_DIR)/, $(chem_tests_:.cc=.o))

