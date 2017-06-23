pssp_tests_basis_=pssp/tests/basis/Test.cc

pssp_tests_basis_SRCS=\
     $(addprefix $(SRC_DIR)/, $(pssp_tests_basis_))
pssp_tests_basis_OBJS=\
     $(addprefix $(BLD_DIR)/, $(pssp_tests_basis_:.cc=.o))

